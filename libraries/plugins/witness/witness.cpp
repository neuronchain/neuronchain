/*
 * Copyright (c) 2015 Cryptonomex, Inc., and contributors.
 *
 * The MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <graphene/witness/witness.hpp>

#include <graphene/chain/database.hpp>
#include <graphene/chain/witness_object.hpp>

#include <graphene/utilities/key_conversion.hpp>

#include <fc/smart_ref_impl.hpp>
#include <fc/thread/thread.hpp>

#include <fc/thread/mutex.hpp>
#include <fc/thread/unique_lock.hpp>

#include <iostream>

using namespace graphene::witness_plugin;
using std::string;
using std::vector;

namespace bpo = boost::program_options;

void new_chain_banner( const graphene::chain::database& db )
{
   std::cerr << "\n"
      "********************************\n"
      "*                              *\n"
      "*   ------- NEW CHAIN ------   *\n"
      "*   - Welcome to Graphene! -   *\n"
      "*   ------------------------   *\n"
      "*                              *\n"
      "********************************\n"
      "\n";
   if( db.get_slot_at_time( fc::time_point::now() ) > 200 )
   {
      std::cerr << "Your genesis seems to have an old timestamp\n"
         "Please consider using the --genesis-timestamp option to give your genesis a recent timestamp\n"
         "\n"
         ;
   }
   return;
}

void witness_plugin::plugin_set_program_options(
   boost::program_options::options_description& command_line_options,
   boost::program_options::options_description& config_file_options)
{
   auto default_priv_key = fc::ecc::private_key::regenerate(fc::sha256::hash(std::string("faucet_account")));
   string witness_id_example = fc::json::to_string(chain::witness_id_type(5));
   command_line_options.add_options()
         ("enable-stale-production", bpo::bool_switch()->notifier([this](bool e){_production_enabled = e;}), "Enable block production, even if the chain is stale.")
         ("required-participation", bpo::bool_switch()->notifier([this](int e){_required_witness_participation = uint32_t(e*GRAPHENE_1_PERCENT);}), "Percent of witnesses (0-99) that must be participating in order to produce blocks")
         ("witness-id,w", bpo::value<vector<string>>()->composing()->multitoken(),
          ("ID of witness controlled by this node (e.g. " + witness_id_example + ", quotes are required, may specify multiple times)").c_str())
         ("private-key", bpo::value<vector<string>>()->composing()->multitoken()->
          DEFAULT_VALUE_VECTOR(std::make_pair(chain::public_key_type(default_priv_key.get_public_key()), graphene::utilities::key_to_wif(default_priv_key))),
          "Tuple of [PublicKey, WIF private key] (may specify multiple times)")
         //("enable-transaction-generator", bpo::bool_switch()->notifier([this](bool e){_generation_enabled = e;}), "Enable transaction generator node.")
         //("generate-transactions", bpo::value<uint32_t>()->default_value(0), "Generate fixed number of transactions periodically")
         ;
   config_file_options.add(command_line_options);
}

std::string witness_plugin::plugin_name()const
{
   return "witness";
}

void witness_plugin::plugin_initialize(const boost::program_options::variables_map& options)
{ try {
   ilog("witness plugin:  plugin_initialize() begin");
   _options = &options;
   LOAD_VALUE_SET(options, "witness-id", _witnesses, chain::witness_id_type)

   if( options.count("private-key") )
   {
      const std::vector<std::string> key_id_to_wif_pair_strings = options["private-key"].as<std::vector<std::string>>();
      for (const std::string& key_id_to_wif_pair_string : key_id_to_wif_pair_strings)
      {
         auto key_id_to_wif_pair = graphene::app::dejsonify<std::pair<chain::public_key_type, std::string> >(key_id_to_wif_pair_string);
         //idump((key_id_to_wif_pair));
         ilog("Public Key: ${public}", ("public", key_id_to_wif_pair.first));
         fc::optional<fc::ecc::private_key> private_key = graphene::utilities::wif_to_key(key_id_to_wif_pair.second);
         if (!private_key)
         {
            // the key isn't in WIF format; see if they are still passing the old native private key format.  This is
            // just here to ease the transition, can be removed soon
            try
            {
               private_key = fc::variant(key_id_to_wif_pair.second).as<fc::ecc::private_key>();
            }
            catch (const fc::exception&)
            {
               FC_THROW("Invalid WIF-format private key ${key_string}", ("key_string", key_id_to_wif_pair.second));
            }
         }
         _private_keys[key_id_to_wif_pair.first] = *private_key;
      }
   }
   ilog("witness plugin:  plugin_initialize() end");
} FC_LOG_AND_RETHROW() }

void witness_plugin::plugin_startup()
{ try {
   ilog("witness plugin:  plugin_startup() begin");
   chain::database& d = database();

   if( !_witnesses.empty() )
   {
      ilog("Launching block production for ${n} witnesses.", ("n", _witnesses.size()));
      app().set_block_production(true);
      if( _production_enabled )
      {
         if( d.head_block_num() == 0 )
            new_chain_banner(d);
         _production_skip_flags |= graphene::chain::database::skip_undo_history_check;
      }
      schedule_production_loop();
   } else
      elog("No witnesses configured! Please add witness IDs and private keys to configuration.");
   //if (_options["generate-transactions"].as<uint32_t>() > 0)
   if (_generation_enabled)
   {
      //schedule_generation_loop();
      fc::time_point now = fc::time_point::now();
      fc::time_point next_wakeup( fc::time_point::now() + fc::seconds( 40 ) );
      //wdump( (now.time_since_epoch().count())(next_wakeup.time_since_epoch().count()) );
      _transaction_generation_task = fc::schedule([this] { transaction_generation_loop(); },
                                                  next_wakeup, "Transaction generator");
   }
   ilog("witness plugin:  plugin_startup() end");
} FC_CAPTURE_AND_RETHROW() }

void witness_plugin::plugin_shutdown()
{
   return;
}

void witness_plugin::schedule_production_loop()
{
   //Schedule for the next second's tick regardless of chain state
   // If we would wait less than 50ms, wait for the whole second.
   fc::time_point now = fc::time_point::now();
   int64_t time_to_next_second = 1000000 - (now.time_since_epoch().count() % 1000000);
   if( time_to_next_second < 50000 )      // we must sleep for at least 50ms
       time_to_next_second += 1000000;

   fc::time_point next_wakeup( now + fc::microseconds( time_to_next_second ) );

   //wdump( (now.time_since_epoch().count())(next_wakeup.time_since_epoch().count()) );
   _block_production_task = fc::schedule([this]{block_production_loop();},
                                         next_wakeup, "Witness Block Production");
}

block_production_condition::block_production_condition_enum witness_plugin::block_production_loop()
{
   try
   {
      chain::database& db = database();
      FC_ASSERT(_witnesses.size() > 0);
      chain::witness_id_type signing_witness = *_witnesses.begin();
      graphene::chain::public_key_type public_key = signing_witness( db ).signing_key;
      auto private_key_itr = _private_keys.find( public_key );
      if( private_key_itr == _private_keys.end() )
      {
         elog("No key for signing packet!");
      }
      auto packet = db.push_packet(signing_witness, private_key_itr->second);
      fc::async( [this,packet](){ p2p_node().broadcast(net::packet_message(packet)); } );
   }
   catch( const fc::canceled_exception& )
   {
      //We're trying to exit. Go ahead and let this one out.
      throw;
   }
   catch( const fc::exception& e )
   {
      elog("Got exception while pushing packet:\n${e}", ("e", e.to_detail_string()));
   }


   block_production_condition::block_production_condition_enum result;
   fc::mutable_variant_object capture;
   try
   {
      result = maybe_produce_block(capture);
   }
   catch( const fc::canceled_exception& )
   {
      //We're trying to exit. Go ahead and let this one out.
      throw;
   }
   catch( const fc::exception& e )
   {
      elog("Got exception while generating block:\n${e}", ("e", e.to_detail_string()));
      result = block_production_condition::exception_producing_block;
   }

   switch( result )
   {
      case block_production_condition::produced: 
         ilog("Generated block #${n} with timestamp ${t} at time ${c}", (capture));
         break;
      
      case block_production_condition::not_synced:
         ilog("Not producing block because production is disabled until we receive a recent block (see: --enable-stale-production)");
         break;
      case block_production_condition::not_my_turn:
         //ilog("Not producing block because it isn't my turn");
         break;
      case block_production_condition::not_time_yet:
         //ilog("Not producing block because slot has not yet arrived");
         break;
      case block_production_condition::no_private_key:
         ilog("Not producing block because I don't have the private key for ${scheduled_key}", (capture) );
         break;
      case block_production_condition::low_participation:
         elog("Not producing block because node appears to be on a minority fork with only ${pct}% witness participation", (capture) );
         break;
      case block_production_condition::lag:
         elog("Not producing block because node didn't wake up within 500ms of the slot time.");
         break;
      case block_production_condition::consecutive:
         elog("Not producing block because the last block was generated by the same witness.\nThis node is probably disconnected from the network so block production has been disabled.\nDisable this check with --allow-consecutive option.");
         break;
      case block_production_condition::exception_producing_block:
         elog( "exception prodcing block" );
         break;
   }

   schedule_production_loop();
   return result;
}

block_production_condition::block_production_condition_enum witness_plugin::maybe_produce_block( fc::mutable_variant_object& capture )
{
   chain::database& db = database();
   fc::time_point now_fine = fc::time_point::now();
   //fc::time_point_sec now = now_fine + fc::microseconds( 500000 );
   fc::time_point_sec now = now_fine + fc::microseconds( 1000000 );

   // If the next block production opportunity is in the present or future, we're synced.
   if( !_production_enabled )
   {
      if( db.get_slot_time(1) >= now )
         _production_enabled = true;
      else
         return block_production_condition::not_synced;
   }

   // is anyone scheduled to produce now or one second in the future?
   uint32_t slot = db.get_slot_at_time( now );
   if( slot == 0 )
   {
      capture("next_time", db.get_slot_time(1));
      return block_production_condition::not_time_yet;
   }

   //
   // this assert should not fail, because now <= db.head_block_time()
   // should have resulted in slot == 0.
   //
   // if this assert triggers, there is a serious bug in get_slot_at_time()
   // which would result in allowing a later block to have a timestamp
   // less than or equal to the previous block
   //
   assert( now > db.head_block_time() );

   graphene::chain::witness_id_type scheduled_witness = db.get_scheduled_witness( slot );
   // we must control the witness scheduled to produce the next block.
   if( _witnesses.find( scheduled_witness ) == _witnesses.end() )
   {
      capture("scheduled_witness", scheduled_witness);
      return block_production_condition::not_my_turn;
   }

   fc::time_point_sec scheduled_time = db.get_slot_time( slot );
   graphene::chain::public_key_type scheduled_key = scheduled_witness( db ).signing_key;
   auto private_key_itr = _private_keys.find( scheduled_key );

   if( private_key_itr == _private_keys.end() )
   {
      capture("scheduled_key", scheduled_key);
      return block_production_condition::no_private_key;
   }

   uint32_t prate = db.witness_participation_rate();
   if( prate < _required_witness_participation )
   {
      capture("pct", uint32_t(100*uint64_t(prate) / GRAPHENE_1_PERCENT));
      return block_production_condition::low_participation;
   }

   if( llabs((scheduled_time - now).count()) > fc::milliseconds( 10000 ).count() )
   {
      capture("scheduled_time", scheduled_time)("now", now);
      return block_production_condition::lag;
   }


   //wlog("Begin block generation");
   //auto start = fc::time_point::now();
   auto block = db.generate_block(
      scheduled_time,
      scheduled_witness,
      private_key_itr->second,
      _production_skip_flags
      );
   //auto dur = fc::time_point::now() - start;
   //idump((dur.count() / 1000));
   ilog("New block with ${txs} transactions", ("txs", block.transactions.size()));
   capture("n", block.block_num())("t", block.timestamp)("c", now);
   //p2p_node().broadcast(net::block_message(block));
   fc::async( [this,block](){ p2p_node().broadcast(net::block_message(block)); } );

  

//genth.async(flood);
   return block_production_condition::produced;
}


void witness_plugin::schedule_generation_loop()
{
   //Schedule for the next second's tick regardless of chain state
   // If we would wait less than 50ms, wait for the whole second.
   static size_t steps = 0;
   fc::time_point now = fc::time_point::now();
   int64_t time_to_next_second = 1000000 - (now.time_since_epoch().count() % 1000000);
   if( time_to_next_second < 50000 )      // we must sleep for at least 50ms
      time_to_next_second += 1000000;

   auto add = fc::seconds(6);
   if (steps++ > 3) {
      steps = 0;
      add = fc::seconds(25);
   }
   fc::time_point next_wakeup( now + fc::microseconds( time_to_next_second) + add );

   //wdump( (now.time_since_epoch().count())(next_wakeup.time_since_epoch().count()) );
   _transaction_generation_task = fc::schedule([this]{transaction_generation_loop();},
                                         next_wakeup, "Transaction generator");
}

void witness_plugin::transaction_generation_loop()
{
   //wlog("transaction generation loop start");
   //static fc::thread genth("transfers_generator");
   using namespace graphene::chain;
   uint32_t skip = graphene::chain::database::skip_witness_signature |
                   graphene::chain::database::skip_transaction_signatures |
                   graphene::chain::database::skip_tapos_check |
                   graphene::chain::database::skip_authority_check |
                   graphene::chain::database::skip_transaction_dupe_check;
   chain::database& db = database();
   //fc::optional<asset_object> asset_obj = asset( 0 );
   static size_t step = 1;
   static uint64_t nonce = 0;
   static uint16_t to = 6;
   //auto flood = [&]{
       //fc::unique_lock<fc::mutex> lock(mx, fc::try_to_lock_t{});
       //fc::unique_lock<fc::mutex> lock(mx, fc::time_point::now() + fc::seconds(0.5));
       //if (!mx.try_lock()) return;
       //fc::time_point begin_gen = fc::time_point::now();
       signed_transaction tx;

       //tx.ref_block_num = db.head_block_num();
       //tx.ref_block_prefix = db.head_block_id();

       account_id_type from_id = account_id_type{17};
       account_id_type to_id = account_id_type{to++};

       transfer_operation xfer_op;

       xfer_op.from = from_id;
       xfer_op.to = to_id;
       xfer_op.amount =asset( 1 );//asset_obj->amount(1);//asset_obj->amount_from_string("1");

       //xfer_op.memo = memo_data();

       //xfer_op.memo->set_message(get_private_key(from_account.options.memo_key), to_account.options.memo_key, memo);

       tx.operations.push_back(xfer_op);
       const auto& fees = *db.get_global_properties().parameters.current_fees;
       for( auto& op : tx.operations )
          fees.set_fee(op);
       //auto& op = xfer_op;
       //fees.set_fee(xfer_op);
       //set_operation_fees( tx, db.get_global_properties().parameters.current_fees);
       //tx.validate();

       auto dynprops = db.get_dynamic_global_properties();
       const auto& chain_id = db.get_chain_id();

       //tx.set_reference_block( dynprops.head_block_id );
       tx.signatures.clear();

       //wdump((_options["generate-transactions"].as<uint32_t>()));
       const int64_t transactions_num = _options["generate-transactions"].as<uint32_t>();//6000;
       std::vector<signed_transaction> trs;
       for (size_t i = 0; i < transactions_num; ++i) {
       //fc::unique_lock<fc::mutex> lock(mx, fc::time_point::now() + fc::seconds(0.1));
          auto exp = /*fc::time_point::now()*/ dynprops.time + fc::seconds(i % 86000 + 10);
          tx.set_expiration(exp);
          tx.operations.front().get<transfer_operation>().amount = asset( (nonce * transactions_num + i) % 100000000 + 1);
          //xfer_op.memo->nonce = nonce++;
          //std::string str = std::to_string(nonce);
          //xfer_op.memo->message = std::vector<char>(str.begin(), str.end());

          //tx = sign_transaction(tx, false);

          //auto res = db.push_transaction( tx, skip );
          auto trx_cpy = tx;
          trx_cpy.sign( _private_keys.begin()->second, chain_id );
          trs.push_back(std::move(trx_cpy));

          //db.bundle_transaction(tx, skip);

          //fc::time_point end_gen = fc::time_point::now();
          //ilog("Generation time: ${gent}", ("gent", end_gen - begin_gen));
          //mx.unlock();
       }
       for (auto& tr : trs)
          db.bundle_transaction(tr);
       ++nonce;
   //};
   //flood();

   if (nonce > 100)
      nonce = 0;
   if (to > 13) to = 6;

   try
   {
      chain::database& db = database();
      //FC_ASSERT(_witnesses.size() > 0);
      chain::witness_id_type signing_witness = *(db.get_global_properties().active_witnesses.begin());//*_witnesses.begin();
      graphene::chain::public_key_type public_key = signing_witness( db ).signing_key;
      auto private_key_itr = _private_keys.find( public_key );
      if( private_key_itr == _private_keys.end() )
      {
         elog("No key for signing packet!");
      }
      auto packet = db.push_packet(signing_witness, private_key_itr->second);
      fc::async( [this,packet](){ p2p_node().broadcast(net::packet_message(packet)); } );
      ilog("Broadcasting packet with ${n} transactions", ("n", packet.transactions.size()));
   }
   catch( const fc::canceled_exception& )
   {
      //We're trying to exit. Go ahead and let this one out.
      throw;
   }
   catch( const fc::exception& e )
   {
      elog("Got exception while pushing packet:\n${e}", ("e", e.to_detail_string()));
   }

   schedule_generation_loop();
}