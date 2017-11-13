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

#include <queue>
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

      int64_t initial_balance = 30000;

      for(size_t i = 0; i < 6; ++i) {
         upgrade_to_lifetime_member(accounts[i]);
         wit_ids[i] = create_witness(accounts[i], pkeys[i]).id;

         transfer(account_id_type(), accounts[i], asset(initial_balance));
         //fund(db.get(accounts[i]), asset(initial_balance));
      }

      BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 30000);

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

      //BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 99994995);

      /*std::cerr << "checkpoint\nscores: ";
      for (const auto& acc : accounts) {
            std::cerr << db.get(db.get(acc).statistics).importance_score << " ";
      }*/

      transfer(account0_id, account5_id, asset(1001));
      generate_blocks(db.get_dynamic_global_properties().next_maintenance_time, true, skip_sigs);
      
      /*std::cerr << "scores: " << db.get(db.get(account_id_type()).statistics).importance_score << " ";
      for (const auto& acc : accounts) {
            std::cerr << db.get(db.get(acc).statistics).importance_score << " ";
      }*/

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

      int64_t initial_balance = 30000;

      for(size_t i = 0; i < 6; ++i) {
      upgrade_to_lifetime_member(accounts[i]);
      //trx.clear();
         wit_ids[i] = create_witness(accounts[i], pkeys[i]).id;

         transfer(account_id_type(), accounts[i], asset(initial_balance));
      }

      BOOST_CHECK(get_balance(account0_id, asset_id_type{}) == 30000);

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

#if 0
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
#endif

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

#if 0
namespace test {
      //using transfer_graph_type = std::unordered_multimap<object_id_type, object_id_type>;
      using transfer_graph_type = std::unordered_multimap<uint64_t, uint64_t>;
      using node_type = transfer_graph_type::key_type;

std::set<node_type> get_gamma(const transfer_graph_type& transfer_graph, const node_type& base_node) 
{
      std::set<node_type> gamma{base_node};
      size_t bucket_num = transfer_graph.bucket(base_node);
      auto it = transfer_graph.cbegin(bucket_num);
      for (; it != transfer_graph.cend(bucket_num); ++it)
            gamma.insert(it->second);
      return gamma;
};

fc::optional<std::unordered_set<node_type>> get_core(const transfer_graph_type& transfer_graph, node_type node, 
                                                      double epsilon, uint32_t mu) 
{
      //size_t bucket_num = transfer_graph.bucket(node);
      //auto it = transfer_graph.cbegin(bucket_num);
    
      //if (transfer_graph.bucket_size(bucket_num) < mu)
      //      return {};

      std::set<node_type> first = get_gamma(transfer_graph, node);
      if (first.size() < mu)
            return {};

      std::unordered_set<node_type> core;
      auto it = first.cbegin();
      //for (; it != transfer_graph.cend(bucket_num); ++it) {
      for (; it != first.cend(); ++it) {
            std::set<node_type> second = get_gamma(transfer_graph, *it);
            std::set<node_type> intersection;
            std::set_intersection(first.cbegin(), first.cend(), second.cbegin(), second.cend(), 
                        std::inserter(intersection, intersection.begin()));
            double similarity = intersection.size() / std::sqrt(first.size() * second.size());

            std::cerr << node << " <-> " << *it << " : " << 
                  intersection.size() << " / √(" << first.size() << " * " << second.size() << ") = " <<
                  similarity << std::endl;

            if (similarity >= epsilon)
                  core.insert(*it);
      }
      if (core.size() >= mu) return core;
      return {};
      //return (sum > mu) ? core : {};
};

std::unordered_map<node_type, uint32_t> cluster_graph_simple(const transfer_graph_type& transfer_graph, double epsilon, uint32_t mu)
{
      std::unordered_map<node_type, uint32_t> node_ids;
      const uint32_t non_member_id = 0;
      const uint32_t outlier_id = 1;
      const uint32_t hub_id = 2;
      uint32_t last_cluster = 2;

      auto skip_same = [](transfer_graph_type::const_iterator it, transfer_graph_type::const_iterator end) {
            transfer_graph_type::const_iterator last = it++;
            for(; it != end && it->first == last->first; last = it++);
            return it;
      };

      auto it = transfer_graph.begin();
      for(; it != transfer_graph.end(); it = skip_same(it, transfer_graph.end())) {
            //std::cerr << it->first.instance() << " ";
            const auto& node = it->first;
            if (node_ids.find(node) != node_ids.end())
                  continue;

            std::queue<node_type> queue;
            auto core_optional = get_core(transfer_graph, node, epsilon, mu);
            if (core_optional.valid()) {
                  last_cluster++;
                  node_ids[node] = last_cluster;
                  for (const auto& neighbour : *core_optional) {
                        /*auto neighbour_id = node_ids.find(neighbour);
                        if (neighbour_id == node_ids.end()) {
                              node_ids[neighbour] = last_cluster;
                        }*/
                        queue.push(neighbour);
                  }
                  while(!queue.empty()) {
                        auto front = queue.front();
                        queue.pop();
                        auto front_core_optional = get_core(transfer_graph, front, epsilon, mu);
                        if (front_core_optional.valid()) {
                              for (const auto& reachable_node : *front_core_optional) {
                                    auto reachable_id = node_ids.find(reachable_node);
                                    if (reachable_id == node_ids.end()) {
                                          node_ids[reachable_node] = last_cluster;
                                          queue.push(reachable_node);
                                    }
                                    else if (reachable_id->second == non_member_id)
                                          reachable_id->second = last_cluster;
                              }
                        }                        
                  }
            }
            else {
                  node_ids[node] = non_member_id;
            }
      } 
      std::cerr << "!\n";
      for (auto& id : node_ids) {
      //it = transfer_graph.begin();
      //for(; it != transfer_graph.end(); it = skip_same(it, transfer_graph.end())) {
            //std::cerr << id.first.instance() << ' ';
            //auto id = std::make_pair(it->first, node_ids.at(it->first));
            if (id.second == non_member_id) {
                  auto gamma = get_gamma(transfer_graph, id.first);
                  /*bool is_hub = std::inner_product(gamma.begin(), gamma.end(), gamma.begin(), false, std::logical_or<bool>{}, 
                        [&](const decltype(*gamma.begin())& lhs, const decltype(*gamma.begin())& rhs) {
                              return node_ids.at(lhs) != node_ids.at(rhs);
                  });*/
                  bool is_hub = false;
                  auto g_it1 = gamma.cbegin();
                  auto g_end = gamma.cend();
                  std::cerr << id.first << " -> " << std::flush;
                  for (const auto& gnode : gamma)
                        std::cerr << gnode << ' ' << std::flush;
                  for (; g_it1 != g_end && !is_hub; ++g_it1) {
                        for (auto g_it2 = std::next(g_it1); g_it2 != g_end && !is_hub; ++g_it2) {
                              auto id1 = node_ids.at(*g_it1);
                              auto id2 = node_ids.at(*g_it2);
                              if (id1 > hub_id && id2 > hub_id && id1 != id2)
                                    is_hub = true;
                        }
                  }
                  std::cerr << "!\n";
                  id.second = is_hub ? hub_id : outlier_id;
            }
      }
      std::cerr << "!\n";
      return node_ids;
}
}
#endif

BOOST_AUTO_TEST_CASE( scan_test )
{
      try {
            using transfer_graph_type = std::unordered_multimap<uint64_t, uint64_t>;
            
            {
                  std::vector<std::vector<uint32_t>> edges{
                        {1,4,5,6},        //0
                        {0,2,5},          //1
                        {1,3,5},          //2
                        {2,4,5,6},        //3
                        {0,3,5,6},        //4
                        {0,1,2,3,4},      //5
                        {0,3,4,7,10,11},  //6
                        {6,8,11,12},      //7
                        {7,9,12},         //8
                        {8,10,12,13},     //9
                        {6,9,11,12},      //10
                        {6,7,10,12},      //11
                        {7,8,9,10,11},    //12
                        {9}               //13
                  };

                  std::vector<uint64_t> accounts(edges.size());
                  std::iota(accounts.begin(), accounts.end(), 0);
                  /*for(size_t i = 0; i < edges.size(); ++i)
                        accounts.emplace_back(account_id_type(i));*/

                  std::unordered_map<uint64_t, uint64_t> accounts_map;
                  for (size_t i = 0; i < accounts.size(); ++i)
                        accounts_map.emplace(i, i);

                  //transfer_graph_type test_graph;
                  boost::numeric::ublas::matrix<double> outlink(edges.size(), edges.size(), 0);
                  double outlink_value = db.get_global_properties().parameters.min_transfer_for_clustering.value + 2;
                  for(size_t i = 0; i < edges.size(); ++i) {
                        auto& edges_row = edges[i];
                        for (auto vertex : edges_row)
                              outlink(accounts[i], accounts[vertex]) = outlink_value;
                              //test_graph.emplace(accounts[i], accounts[vertex]);
                  }

                  /*auto it = test_graph.begin();
                  for (; it != test_graph.end(); ++it)
                        std::cerr << it->first << " -> " << it->second << std::endl;
                  std::cerr << '\n';*/

                  auto min_transfer = db.get_global_properties().parameters.min_transfer_for_clustering.value;
                  auto clusters = graphene::chain::detail::cluster_graph_simple(outlink, accounts_map, 0.7, 3, min_transfer);
                  //auto clusters = test::cluster_graph_simple(test_graph, 0.7, 3);
            
                  for (const auto& p : clusters)
                        std::cerr << p.first << " -> " << p.second << std::endl;

                  std::vector<uint32_t> cluster1{0,1,2,3,4,5};
                  std::vector<uint32_t> cluster2{7,8,9,10,11,12};

                  for (auto k : cluster1)
                        BOOST_CHECK(clusters[accounts[cluster1[0]]] == clusters[accounts[k]]);
                  for (auto k : cluster2)
                        BOOST_CHECK(clusters[accounts[cluster2[0]]] == clusters[accounts[k]]);
                  BOOST_CHECK(clusters[accounts[cluster1[0]]] != clusters[accounts[cluster2[0]]]);

                  BOOST_CHECK(clusters[accounts[6]] == 2);  //hub
                  BOOST_CHECK(clusters[accounts[13]] == 1); //outlier     
            }

            {
                  std::vector<std::vector<uint32_t>> edges{
                        {8},                    //0
                        {2,4,6,8},              //1
                        {1,4,6,7,8,13},         //2
                        {4,7,8,13,14},          //3
                        {1,2,3,5,6,7,8,13,14,15},//4
                        {4,15,16},              //5
                        {1,2,4,7,13,14,15},       //6
                        {2,3,4,6,8,15},         //7
                        {0,1,2,3,4,7},          //8
                        {10,11,12},             //9
                        {9,11,12,16},           //10
                        {9,10,12,16},           //11
                        {9,10,11},              //12
                        {2,3,4,6,14},           //13
                        {3,4,6,13},             //14
                        {4,5,6,7},              //15
                        {5,10,11}               //16
                  };
                  /*std::vector<std::vector<uint32_t>> edges{
                        {8},                    //0
                        {2,4,6,8},              //1
                        {4,6,7,8,13},         //2
                        {4,7,8,13,14},          //3
                        {5,6,7,13,14,15},     //4
                        {15,16},              //5
                        {7,13,14,15},       //6
                        {8,15},         //7
                        {},          //8
                        {10,11,12},             //9
                        {11,12},              //10
                        {12},              //11
                        {},              //12
                        {14},           //13
                        {},             //14
                        {},              //15
                        {}               //16
                  };*/

                  std::vector<uint64_t> accounts(edges.size());
                  std::iota(accounts.begin(), accounts.end(), 0);
                  /*for(size_t i = 0; i < edges.size(); ++i)
                        accounts.emplace_back(account_id_type(i));*/

                  std::unordered_map<uint64_t, uint64_t> accounts_map;
                  for (size_t i = 0; i < accounts.size(); ++i)
                        accounts_map.emplace(i, i);
                  //transfer_graph_type test_graph;
                  boost::numeric::ublas::matrix<double> outlink(edges.size(), edges.size(), 0);
                  auto outlink_value = db.get_global_properties().parameters.min_transfer_for_clustering.value + 1;
                  for(size_t i = 0; i < edges.size(); ++i) {
                        auto& edges_row = edges[i];
                        for (auto vertex : edges_row) {
                              /*outlink(accounts[i], accounts[vertex]) = outlink_value;
                              test_graph.emplace(accounts[i], accounts[vertex]);*/
                              outlink(i, vertex) = outlink_value;
                              //test_graph.emplace(i, vertex);
                        }
                  }

                  /*std::cerr << '\n';
                  auto it = test_graph.begin();
                  for (; it != test_graph.end(); ++it)
                        std::cerr << it->first << " -> " << it->second << std::endl;
                  std::cerr << '\n';*/

                  auto min_transfer = db.get_global_properties().parameters.min_transfer_for_clustering.value;
                  auto clusters = graphene::chain::detail::cluster_graph_simple(outlink, accounts_map,  0.671, 2, min_transfer);
                  //auto clusters2 = test::cluster_graph_simple(test_graph, 0.671, 2);
            
                  for (const auto& p : clusters)
                        std::cerr << p.first << " -> " << p.second << std::endl;

                  std::vector<uint32_t> cluster1{1,2,3,4,6,7,8,13,14,15};
                  std::vector<uint32_t> cluster2{9,10,11,12};

                  for (auto k : cluster1)
                        BOOST_CHECK(clusters[accounts[cluster1[0]]] == clusters[accounts[k]]);
                  for (auto k : cluster2)
                        BOOST_CHECK(clusters[accounts[cluster2[0]]] == clusters[accounts[k]]);
                  BOOST_CHECK(clusters[accounts[cluster1[0]]] != clusters[accounts[cluster2[0]]]);

                  BOOST_CHECK(clusters[accounts[5]] != clusters[accounts[cluster1[0]]]);
                  //BOOST_CHECK(clusters[accounts[5]] != clusters[accounts[16]]);

                  BOOST_CHECK(clusters[accounts[16]] != clusters[accounts[cluster2[0]]]);

                  BOOST_CHECK(clusters[accounts[0]] == 1); //outlier
            }
       

      } FC_LOG_AND_RETHROW()
}


BOOST_AUTO_TEST_CASE( rank_test )
{
      try {
            {
                  std::vector<std::vector<uint32_t>> edges{
                        {1},                    //0
                        {2},                    //1
                        {0,3,6},                //2
                        {4},                    //3
                        {5},                    //4
                        {3},                    //5
                        {},                     //6
                  };

                  std::vector<uint64_t> accounts(edges.size());
                  std::iota(accounts.begin(), accounts.end(), 0);
                  /*for(size_t i = 0; i < edges.size(); ++i)
                        accounts.emplace_back(account_id_type(i));*/

                  std::unordered_map<uint64_t, uint64_t> accounts_map;
                  for (size_t i = 0; i < accounts.size(); ++i)
                        accounts_map.emplace(i, i);

                  //transfer_graph_type test_graph;
                  boost::numeric::ublas::matrix<double> outlink(edges.size(), edges.size());
                  auto outlink_value = db.get_global_properties().parameters.min_transfer_for_clustering.value + 1;
                  for(size_t i = 0; i < edges.size(); ++i) {
                        auto& edges_row = edges[i];
                        for (auto vertex : edges_row)
                              outlink(accounts[i], accounts[vertex]) = outlink_value;
                              //test_graph.emplace(accounts[i], accounts[vertex]);
                  }
                  //for (size_t j = 0; j < outlink.size2(); ++j) {
                  for (auto it = outlink.begin1(); it != outlink.end1(); ++it) {
                        //boost::numeric::ublas::matrix_column<decltype(outlink)> col(outlink, j);
                        double norm = std::accumulate(it.begin(), it.end(), 0.);
                        if (norm > 0)
                           for (auto it2 = it.begin(); it2 != it.end(); ++it2)
                                 *it2 /= norm;
                  }

                  /*auto it = test_graph.begin();
                  for (; it != test_graph.end(); ++it)
                        std::cerr << it->first.instance() << " -> " << it->second.instance() << std::endl;
                  std::cerr << '\n';*/

                  auto min_transfer = db.get_global_properties().parameters.min_transfer_for_clustering.value;
                  double mu = db.get_global_properties().parameters.clustering_mu;
                  double epsilon = db.get_global_properties().parameters.clustering_epsilon;
                  auto clusters = graphene::chain::detail::cluster_graph_simple(outlink, accounts_map, 0.7, 2, min_transfer);

                  double r_mu = db.get_global_properties().parameters.rank_mu;
                  double r_etha = db.get_global_properties().parameters.rank_etha;
                  double r_epsilon = db.get_global_properties().parameters.rank_epsilon;

                  graphene::chain::detail::getNCDAwareRank(outlink, clusters, r_etha, r_mu, r_epsilon);                  
            }


      } FC_LOG_AND_RETHROW()
}
BOOST_AUTO_TEST_SUITE_END()   
