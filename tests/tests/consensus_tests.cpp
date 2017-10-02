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

#include <boost/test/unit_test.hpp>

#include <graphene/chain/database.hpp>
#include <graphene/chain/account_object.hpp>
#include <graphene/chain/witness_object.hpp>

#include <graphene/utilities/tempdir.hpp>   

#include <fc/crypto/digest.hpp>

#include <random>

#include "../common/database_fixture.hpp"

using namespace graphene::chain;
using namespace graphene::chain::test;

BOOST_FIXTURE_TEST_SUITE( consensus_tests, database_fixture )

BOOST_AUTO_TEST_CASE( consensus_test )
{
   try {
      ACTORS( (account0)(account1)(account2)(account3)(account4)(account5) );

      generate_block();

      account_id_type accounts[] = {account0_id, account1_id, account2_id,
                                       account3_id, account4_id, account5_id};

      fc::ecc::private_key pkeys[] = {account0_private_key, account1_private_key, account2_private_key,
					account3_private_key, account4_private_key, account5_private_key};

      witness_id_type wit_ids[6];

      int64_t initial_balance = 100000000;

      for(size_t i = 0; i < 6; ++i) {
         upgrade_to_lifetime_member(accounts[i]);
         wit_ids[i] = create_witness(accounts[i], pkeys[i]).id;

         transfer(account_id_type(), accounts[i], asset(initial_balance));
      }

      BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 100000000);

      for (size_t i = 0; i < 5; ++i)
         for (size_t j = i + 1; j < 6; ++j)
            transfer(accounts[i], accounts[j], asset(1001));

      db.modify( db.get_global_properties(), [&]( global_property_object& _gpo )
      {
         _gpo.parameters.maximum_witness_count = 1;
      } );

      auto skip_sigs = database::skip_transaction_signatures | database::skip_authority_check; 

      //generate_blocks(db.get_dynamic_global_properties().next_maintenance_time + fc::seconds(5), true, skip_sigs);
      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);

      BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 99994995);

      transfer(account0_id, account5_id, asset(1001));
      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);

      const auto& witnesses = db.get_global_properties().active_witnesses;

      // make sure we're in active_witnesses
      auto itr = std::find(witnesses.begin(), witnesses.end(), wit_ids[0]);
      BOOST_CHECK(itr != witnesses.end());

      itr = std::find(witnesses.begin(), witnesses.end(), wit_ids[1]);
      BOOST_CHECK(itr == witnesses.end());

      auto& s0 = account0_id(db).statistics(db);
      //std::cout << account0_id(db).name << " -> " << uint64_t(account0_id) << " :: " << s0.importance_score << '\n';
      BOOST_CHECK(s0.importance_score == 6);

      } FC_LOG_AND_RETHROW()
}

BOOST_AUTO_TEST_CASE( cycles_test )
{
   try {
      ACTORS( (account0)(account1)(account2)(account3)(account4)(account5) );

      account_id_type accounts[] = {account0_id, account1_id, account2_id,
                                       account3_id, account4_id, account5_id};

      fc::ecc::private_key pkeys[] = {account0_private_key, account1_private_key, account2_private_key,
               account3_private_key, account4_private_key, account5_private_key};

      witness_id_type wit_ids[6];

      int64_t initial_balance = 100000000;

      for(size_t i = 0; i < 6; ++i) {
      upgrade_to_lifetime_member(accounts[i]);
      //trx.clear();
         wit_ids[i] = create_witness(accounts[i], pkeys[i]).id;

         transfer(account_id_type(), accounts[i], asset(initial_balance));
      }

      BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 100000000);

      for (size_t i = 0; i < 5; ++i)
         for (size_t j = i + 1; j < 6; ++j)
            transfer(accounts[i], accounts[j], asset(1001));

      transfer(account2_id, account0_id, asset(1001));
      //transfer(account5_id, account0_id, asset(1001));

      db.modify( db.get_global_properties(), [&]( global_property_object& _gpo )
      {
         _gpo.parameters.maximum_witness_count = 1;
      } );

      auto skip_sigs = database::skip_transaction_signatures | database::skip_authority_check;
      //BOOST_TEST_MESSAGE("Wait for maintenance interval");
      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, false, skip_sigs);
      //const auto& witnesses = db.get_global_properties().active_witnesses;

      auto& s0 = account0_id(db).statistics(db);
      //std::cout << account0_id(db).name << " -> " << uint64_t(account0_id) << " :: " << s0.importance_score << '\n';
      BOOST_CHECK(s0.importance_score == 3);

      } FC_LOG_AND_RETHROW() 
}

BOOST_AUTO_TEST_CASE( random_test )
{
   try {
      auto skip_sigs = database::skip_transaction_signatures | database::skip_authority_check; 

      std::vector<account_id_type> accounts;
      std::vector<fc::ecc::private_key> pkeys;
      std::vector<witness_id_type> wits;

      
      std::cout << "*" << std::endl;
      int64_t initial_balance = 10000000;
      const size_t accounts_number = 15;
      for (size_t i = 0; i < accounts_number; ++i) {
         std::string name = "accounts" + std::to_string(i);
         pkeys.push_back(generate_private_key(name));
         public_key_type public_key = pkeys[i].get_public_key();
         auto& acc = create_account(name, public_key);
         accounts.push_back(acc.id);
         trx.clear();

         upgrade_to_lifetime_member(accounts[i]);
         wits.push_back(create_witness(accounts[i], pkeys[i]).id);
         trx.clear();

         transfer(account_id_type(), accounts[i], asset(initial_balance));
         trx.clear();

         {
            account_update_operation op;
            op.account = accounts[i];
            op.new_options = accounts[i](db).options;
            op.new_options->votes.insert(wits[i](db).vote_id);// = {wits[i](db).vote_id};
            op.new_options->num_witness = 1;
            op.new_options->num_committee = std::count_if(op.new_options->votes.begin(), op.new_options->votes.end(),
                                                          [](vote_id_type id) { return id.type() == vote_id_type::committee; });
            trx.operations.push_back(op);
            sign( trx, pkeys[i] );
            PUSH_TX( db, trx );
            trx.clear();
         }

         //generate_block(skip_sigs);
      }

      std::cout << "**" << std::endl;
      //BOOST_TEST_MESSAGE("Accounts created");

      generate_block(skip_sigs);

      BOOST_CHECK(get_balance(accounts[0], asset_id_type{}) == 10000000);

      std::random_device r;
      std::default_random_engine e1(r());
      std::uniform_int_distribution<int> uniform_dist(0, accounts_number - 1);

      std::unordered_map<object_id_type, uint32_t> account_transfers;

      std::cout << "***" << std::endl;
      for (size_t i = 0; i < 30; ++i) {
         int from = uniform_dist(e1);
         int to   = uniform_dist(e1);
         if (from == to) continue;

         // std::cout << from << " -> " << to << '\n';
         transfer(accounts[from], accounts[to], asset(1001));
         ++account_transfers[object_id_type(accounts[from])];
         //trx.clear();

         //if (i % 100) generate_block();
         //trx.clear();
//         if (i % 50)
//            generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);
      }
      std::cout << "****" << std::endl;

      //generate_blocks(db.get_dynamic_global_properties().next_maintenance_time + fc::seconds(3), true, skip_sigs);
      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);
      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);

      std::cout << "*****" << std::endl;

      const auto& witnesses = db.get_global_properties().active_witnesses;

      for (const auto& p : account_transfers) {
         std::cout << uint64_t(p.first) << " -> " << p.second << '\n';
      }

      for (const auto& wit : witnesses) {
         auto it = std::find(wits.begin(), wits.end(), wit);
         if (it == wits.end()) std::cout << "Error" << std::endl;
         uint num = it - wits.begin();
         std::cout << "#" << num  << ": " << uint64_t(*it) << " -> " << accounts[num](db).statistics(db).importance_score << " ~ " << account_transfers[object_id_type(accounts[num])] << '\n';
      }

      //auto& s0 = account0_id(db).statistics(db);
      //std::cout << account0_id(db).name << " -> " << uint64_t(account0_id) << " :: " << s0.importance_score << '\n';
      //BOOST_CHECK(s0.importance_score == 6);

      } FC_LOG_AND_RETHROW()
}

#if 0
genesis_state_type make_genesis() {
   genesis_state_type genesis_state;

   genesis_state.initial_timestamp = time_point_sec( GRAPHENE_TESTING_GENESIS_TIMESTAMP );

   auto init_account_priv_key = fc::ecc::private_key::regenerate(fc::sha256::hash(string("null_key")));
   genesis_state.initial_active_witnesses = 1;
   for( int i = 0; i < genesis_state.initial_active_witnesses; ++i )
   {
      auto name = "init"+fc::to_string(i);
      genesis_state.initial_accounts.emplace_back(name,
                                                  init_account_priv_key.get_public_key(),
                                                  init_account_priv_key.get_public_key(),
                                                  true);
      genesis_state.initial_committee_candidates.push_back({name});
      genesis_state.initial_witness_candidates.push_back({name, init_account_priv_key.get_public_key()});
   }
   genesis_state.initial_parameters.current_fees->zero_all_fees();
   genesis_state.immutable_parameters.min_witness_count = 1;
   return genesis_state;
}    

BOOST_AUTO_TEST_CASE( vote_test )
{
   try {
      fc::temp_directory data_dir( graphene::utilities::temp_directory_path() ); 
      database db;
      db.open(data_dir.path(), make_genesis ); 

      ACTORS( (account0)(account1)(account2)(account3)(account4)(account5) );

      generate_block();

      account_id_type accounts[] = {account0_id, account1_id, account2_id,
                                       account3_id, account4_id, account5_id};

      fc::ecc::private_key pkeys[] = {account0_private_key, account1_private_key, account2_private_key,
               account3_private_key, account4_private_key, account5_private_key};

      witness_id_type wit_ids[6];

      int64_t initial_balance = 100000000;

      for(size_t i = 0; i < 6; ++i) {
         upgrade_to_lifetime_member(accounts[i]);
         wit_ids[i] = create_witness(accounts[i], pkeys[i]).id;

         transfer(account_id_type(), accounts[i], asset(initial_balance));

         {
            account_update_operation op;
            op.account = accounts[i];
            op.new_options = accounts[i](db).options;
            op.new_options->votes = {wit_ids[i](db).vote_id};//.insert(wit_ids[1](db).vote_id);
            op.new_options->num_witness = 1;
            op.new_options->num_committee = std::count_if(op.new_options->votes.begin(), op.new_options->votes.end(),
                                                          [](vote_id_type id) { return id.type() == vote_id_type::committee; });
            trx.operations.push_back(op);
            sign( trx, pkeys[i] );
            PUSH_TX( db, trx );
            trx.clear();
         }
      }

      BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 100000000);

      for (size_t i = 1; i < 5; ++i)
         for (size_t j = i + 1; j < 6; ++j)
            transfer(accounts[i], accounts[j], asset(1001));

      db.modify( db.get_global_properties(), [&]( global_property_object& _gpo )
      {
         _gpo.parameters.maximum_witness_count = 1;
      } );

      auto skip_sigs = database::skip_transaction_signatures | database::skip_authority_check; 


     /* {
         account_update_operation op;
         op.account = account1_id;
         op.new_options = account1_id(db).options;
         op.new_options->votes = {wit_ids[1](db).vote_id};//.insert(wit_ids[1](db).vote_id);
         op.new_options->num_witness = 1;
         op.new_options->num_committee = std::count_if(op.new_options->votes.begin(), op.new_options->votes.end(),
                                                       [](vote_id_type id) { return id.type() == vote_id_type::committee; });
         trx.operations.push_back(op);
         sign( trx, pkeys[1] );
         PUSH_TX( db, trx );
         trx.clear();
      }

      {
         account_update_operation op;
         op.account = account2_id;
         op.new_options = account2_id(db).options;
         op.new_options->votes = {wit_ids[2](db).vote_id};//.insert(wit_ids[2](db).vote_id);
         op.new_options->num_witness = 1;
         op.new_options->num_committee = std::count_if(op.new_options->votes.begin(), op.new_options->votes.end(),
                                                       [](vote_id_type id) { return id.type() == vote_id_type::committee; });
         trx.operations.push_back(op);
         sign( trx, pkeys[2] );
         PUSH_TX( db, trx );
         trx.clear();
      }*/

      transfer(account1_id, account5_id, asset(1001));
      transfer(account_id_type(), account1_id, asset(2002));

      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);

      BOOST_CHECK(get_balance(account1_id, asset_id_type{}) == get_balance(account2_id, asset_id_type{}));

      const auto& witnesses = db.get_global_properties().active_witnesses;

      std::cout << witnesses.size() << "!!!\n";

      // make sure we're in active_witnesses
      auto itr = std::find(witnesses.begin(), witnesses.end(), wit_ids[1]);
      BOOST_CHECK(itr != witnesses.end());

      itr = std::find(witnesses.begin(), witnesses.end(), wit_ids[2]);
      BOOST_CHECK(itr == witnesses.end());

      //std::cout << account0_id(db).name << " -> " << uint64_t(account0_id) << " :: " << s0.importance_score << '\n';

      } FC_LOG_AND_RETHROW()
}
#endif

BOOST_AUTO_TEST_SUITE_END()   
