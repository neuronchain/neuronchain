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
#include <graphene/chain/protocol/protocol.hpp>

#include <graphene/chain/account_object.hpp>
#include <graphene/chain/asset_object.hpp>
#include <graphene/chain/proposal_object.hpp>

#include <graphene/db/simple_index.hpp>

#include <fc/crypto/digest.hpp>
#include "../common/database_fixture.hpp"

using namespace graphene::chain;

BOOST_FIXTURE_TEST_SUITE( performance_tests, database_fixture )

BOOST_AUTO_TEST_CASE( sigcheck_benchmark )
{
   fc::ecc::private_key nathan_key = fc::ecc::private_key::generate();
   auto digest = fc::sha256::hash("hello");
   auto sig = nathan_key.sign_compact( digest );
   auto start = fc::time_point::now();
   for( uint32_t i = 0; i < 100000; ++i )
      auto pub = fc::ecc::public_key( sig, digest );
   auto end = fc::time_point::now();
   auto elapsed = end-start;
   wdump(((100000.0*1000000.0) / elapsed.count()) );
}

BOOST_AUTO_TEST_CASE( trs_benchmark )
try {
   fc::ecc::private_key nathan_key = fc::ecc::private_key::generate();
   //const auto& key = register_key(nathan_key.get_public_key());
   private_key_type sam_key = generate_private_key("sam");
   const auto& committee_account = account_id_type()(db);   
   std::vector<account_object> accs;
   for( uint32_t i = 0; i < 1000; ++i )
      accs.emplace_back(create_account("a"+fc::to_string(i), sam_key));
   auto start = fc::time_point::now();
   for( uint32_t i = 0; i < 1e6; ++i )
   {
      /*uint32_t num = i % 990;
      std::vector<account_id_type> tos(10);
      for (uint32_t j = 0; j < 10; ++j)
          tos[j] = accs[num + j].id;*/
      const auto& a = accs[i % 1000];
      fc::async([&]{ 
          transfer( committee_account, a, asset(1 + i) );
      });
      //transfer_to_multiple( committee_account.id, std::move(tos), asset(1000))
      if (i % 100000)
        generate_block();
   }
   
   auto end = fc::time_point::now();
   auto elapsed = end - start;
   wdump( (elapsed) );
   wdump( ((1e6*1e6) / elapsed.count()) );
} FC_LOG_AND_RETHROW()

BOOST_AUTO_TEST_CASE( transfer_benchmark )
{
   fc::ecc::private_key nathan_key = fc::ecc::private_key::generate();
   //const auto& key = register_key(nathan_key.get_public_key());
   private_key_type sam_key = generate_private_key("sam");
   const auto& committee_account = account_id_type()(db);
   auto start = fc::time_point::now();
   for( uint32_t i = 0; i < 1000*1000; ++i )
   {
      const auto& a = create_account("a"+fc::to_string(i), sam_key);
      transfer( committee_account, a, asset(1000) );
   }
   auto end = fc::time_point::now();
   auto elapsed = end - start;
   wdump( (elapsed) );
}


BOOST_AUTO_TEST_SUITE_END()

//#define BOOST_TEST_MODULE "C++ Unit Tests for Graphene Blockchain Database"
#include <cstdlib>
#include <iostream>
#include <boost/test/included/unit_test.hpp>

boost::unit_test::test_suite* init_unit_test_suite(int argc, char* argv[]) {
   std::srand(time(NULL));
   std::cout << "Random number generator seeded to " << time(NULL) << std::endl;
   return nullptr;
}
