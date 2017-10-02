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

#include <fc/crypto/digest.hpp>

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
      const auto& witnesses = db.get_global_properties().active_witnesses;

      auto& s0 = account0_id(db).statistics(db);
      //std::cout << account0_id(db).name << " -> " << uint64_t(account0_id) << " :: " << s0.importance_score << '\n';
      BOOST_CHECK(s0.importance_score == 3);

      // make sure we're not in active_witnesses
      //auto itr = std::find(witnesses.begin(), witnesses.end(), wit_ids[0]);
      //BOOST_CHECK(itr == witnesses.end());

      //BOOST_CHECK(account0.statistics(db).transfers_chronology.rbegin()->second == 0);
      } FC_LOG_AND_RETHROW() 
}

BOOST_AUTO_TEST_SUITE_END()   
