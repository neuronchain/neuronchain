#include <queue>

#include <boost/multiprecision/integer.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/functional/hash.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <fc/smart_ref_impl.hpp>
#include <fc/uint128.hpp>

#include <graphene/chain/database.hpp>

namespace graphene { namespace chain {

using transfer_graph_type = std::unordered_multimap<uint64_t, uint64_t>;
using node_type = transfer_graph_type::key_type;

using matrix_type = boost::numeric::ublas::matrix<double>;
using outlink_matrix_type = matrix_type;
//using outlink_matrix_type = boost::numeric::ublas::compressed_matrix<double>;
//using outlink_matrix_type = std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>>;

/*struct get_transfer_accounts_visitor
{
   struct tranfser_accounts_type
   {
      account_id_type from;
      account_id_type to;
      asset amount;
   };

   typedef fc::optional<tranfser_accounts_type> result_type;

   result_type operator()( const transfer_operation& op )
   {

      return result_type{tranfser_accounts_type{op.from, op.to, op.amount}};
   }

   template <typename Operation>
   result_type operator()( const Operation& op)
   {
      return result_type{};
   }
};*/

outlink_matrix_type database::construct_outlink_matrix(fc::time_point_sec last_maintenance, std::unordered_map<uint64_t, uint64_t>& accounts_map)//, std::unordered_map<object_id_type, uint32_t>& account_transfers)
{
   /*auto transaction_expiration = get_global_properties().parameters.maximum_time_until_expiration;
   if (last_maintenance > head_block_time() - transaction_expiration) {
      transfer_graph_type transfer_graph;

      const auto& all_transactions = get_index_type<transaction_index>().indices().get<by_trx_id>();
      for (const auto& tr : all_transactions) {
         for ( const auto& op : tr.trx.operations) {
            get_transfer_accounts_visitor vtor;
            auto account_pair = op.visit(vtor);
            if (account_pair.valid()) {
               account_id_type from {account_pair->from};
               account_id_type to   {account_pair->to};

               //std::cout << from(*this).name << " : " << uint64_t(from) << " -> " << to(*this).name << '\n';

               auto pair = std::pair<const object_id_type, object_id_type>(to, from); //not a bug
               if (std::find(transfer_graph.cbegin(), transfer_graph.cend(), pair) == transfer_graph.cend())
                  transfer_graph.insert(std::move(pair));

               ++account_transfers[object_id_type{from}];
            }
         }
      }
      return transfer_graph;
   }*/

   auto optional_block = fetch_block_by_id(head_block_id());

   std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> weight_matrix;
   //std::unordered_map<uint64_t, std::unordered_map<uint64_t, double>> outlink_matrix;

   std::unordered_set<uint64_t> accounts_set;

   size_t blocks = 0;
   auto max_blocks = get_global_properties().parameters.importance_score_block_count;
   asset min_transfer_amount = get_global_properties().parameters.min_transfer_for_importance;
   auto min_balance_amount = get_global_properties().parameters.min_balance_for_importance;

    const double log09 = std::log(0.9);
    const uint32_t blocks_per_day = 60 * 60 * 24 / get_global_properties().parameters.block_interval;

    uint32_t head_height = head_block_num();
    uint32_t block_height = head_height;

   while(optional_block.valid() && (blocks++ < max_blocks)) {
      auto& block = *optional_block;    

      for (const auto& tr : block.transactions) {
         for ( const auto& op : tr.operations) {
            get_transfer_accounts_visitor vtor;
            auto transfer_object = op.visit(vtor);
            if (transfer_object.valid() && transfer_object->amount > min_transfer_amount) {
               uint64_t from {object_id_type(transfer_object->from).instance()};
               uint64_t to   {object_id_type(transfer_object->to).instance()};

               auto check_balance = [&](uint64_t acc) -> bool {
                  account_id_type acc_id(acc);
                  auto current_balance = get_balance(acc_id, asset_id_type()).amount.value;
                  return current_balance > min_balance_amount;
               };

               if (check_balance(from) && check_balance(to)) {
                  int64_t floor = std::floor((head_height - block_height) * 1. / blocks_per_day);
                  double weight = transfer_object->amount.amount.value * std::exp(log09 * floor);

                  //std::cout << from(*this).name << " : " << uint64_t(from) << " -> " << to(*this).name << '\n';
                  //auto pair = std::pair<const uint64_t, uint64_t>(from.instance.value, to.instance.instance.value);
                  
                  weight_matrix[to][from] += weight;

                  accounts_set.insert(from);
                  accounts_set.insert(to);
               }
            }
         }
      }

     optional_block = fetch_block_by_id(block.previous);
     block_height = graphene::chain::block_header::num_from_id(block.previous);
   }

   std::vector<uint64_t> accounts(accounts_set.begin(), accounts_set.end());
   std::sort(accounts.begin(), accounts.end());
   // map from account number to continuous index
   //std::unordered_map<uint64_t, uint64_t> accounts_map;
   for (size_t i = 0; i < accounts.size(); ++i)
      accounts_map[accounts.at(i)] = i;

   outlink_matrix_type result(accounts.size(), accounts.size(), 0.);
   for (auto it = weight_matrix.begin(), end = weight_matrix.end(); it != end; ++it) {
      
      size_t row_size = it->second.size();

      for (const auto& tx : it->second) {
            uint64_t to = it->first, from = tx.first;            
            double first_weight = tx.second, second_weight = 0;
            auto from_it = weight_matrix.find(from);
            if (from_it != weight_matrix.end()) {
                  auto to_it = from_it->second.find(to);
                  if (to_it != from_it->second.end())
                        second_weight = to_it->second;
            }
      /*for (const auto& tx : it->second) {
            uint64_t from = it->first, to = tx.first;            
            double first_weight = tx.second, second_weight = 0;
            auto to_it = weight_matrix.find(to);
            if (to_it != weight_matrix.end()) {
                  auto from_it = to_it->second.find(from);
                  if (from_it != to_it->second.end())
                        second_weight = from_it->second;
            }*/

            double outlink = first_weight - second_weight;
            //double outlink = second_weight - first_weight;
            //outlink_matrix[to][from] = outlink > 0 ? outlink : 0;
            result.insert_element(accounts_map.at(from), accounts_map.at(to), outlink > 0 ? outlink : 0);
      }
   }

   //std::cerr << "Res: " << result << "\n!!!\n";

   /*using boost::numeric::ublas::norm_1;
   for (size_t j = 0; j < result.size2(); ++j) {
         boost::numeric::ublas::matrix_column<outlink_matrix_type> col(result, j);
         double norm = norm_1(col);
         if (norm > 0)
            for (size_t i = 0; i < col.size(); ++i)
                  col(i) = col(i) / norm;
   }*/

   return result;
}

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

            /*std::cerr << node << " <-> " << *it << " : " << 
                  intersection.size() << " / √(" << first.size() << " * " << second.size() << ") = " <<
                  similarity << std::endl;*/

            if (similarity >= epsilon)
                  core.insert(*it);
      }
      if (core.size() >= mu) return core;
      return {};
      //return (sum > mu) ? core : {};
};

#if 0
std::unordered_map<node_type, uint32_t> database::detail::cluster_graph_simple(const transfer_graph_type& transfer_graph, double epsilon, uint32_t mu)
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
#endif

const uint32_t MIN_CLUSTER = 2;

std::unordered_map<node_type, uint32_t> detail::cluster_graph_simple(const outlink_matrix_type& outlink_matrix, const std::unordered_map<uint64_t, uint64_t>& accounts_map, double epsilon, uint32_t mu, share_type min_clustering_transfer)
{
      std::unordered_map<node_type, uint32_t> node_ids;
      //const uint32_t non_member_id = 0;
      const uint32_t non_member_id = std::numeric_limits<uint32_t>::max();
      const uint32_t outlier_id = 0;
      const uint32_t hub_id = 1;
      uint32_t last_cluster = 1;

      //auto min_clustering_transfer = get_global_properties().parameters.min_transfer_for_clustering;
      transfer_graph_type transfer_graph;
      for (size_t i = 0; i < outlink_matrix.size1(); ++i) {
            for (size_t j = i + 1; j < outlink_matrix.size2(); ++j)
                  if (outlink_matrix(i,j) + outlink_matrix(j,i) > 0) {//min_clustering_transfer) {
                        transfer_graph.emplace(i, j);
                        transfer_graph.emplace(j, i);
                  }
      }

      auto skip_same = [](transfer_graph_type::const_iterator it, transfer_graph_type::const_iterator end) {
            transfer_graph_type::const_iterator last = it++;
            for(; it != end && it->first == last->first; last = it++);
            return it;
      };

      for(const auto& p: accounts_map) {
            const auto node = p.second;
            if (node_ids.find(node) != node_ids.end())
                  continue;

            std::queue<node_type> queue;
            auto core_optional = get_core(transfer_graph, node, epsilon, mu);
            if (core_optional.valid()) {
                  last_cluster++;
                  node_ids[node] = last_cluster;
                  for (const auto& neighbour : *core_optional) {
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
      for (auto& id : node_ids) {
            if (id.second == non_member_id) {
                  auto gamma = get_gamma(transfer_graph, id.first);
                  /*bool is_hub = std::inner_product(gamma.begin(), gamma.end(), gamma.begin(), false, std::logical_or<bool>{}, 
                        [&](const decltype(*gamma.begin())& lhs, const decltype(*gamma.begin())& rhs) {
                              return node_ids.at(lhs) != node_ids.at(rhs);
                  });*/
                  bool is_hub = false;
                  auto g_it1 = gamma.cbegin();
                  auto g_end = gamma.cend();
                  for (; g_it1 != g_end && !is_hub; ++g_it1) {
                        for (auto g_it2 = std::next(g_it1); g_it2 != g_end && !is_hub; ++g_it2) {
                              auto id1 = node_ids.at(*g_it1);
                              auto id2 = node_ids.at(*g_it2);
                              if (id1 > hub_id && id2 > hub_id && id1 != id2)
                                    is_hub = true;
                        }
                  }
                  id.second = is_hub ? hub_id : outlier_id;
            }
      }
      return node_ids;
}


boost::numeric::ublas::vector<double> detail::getNCDAwareRank(const outlink_matrix_type& outlink_matrix, const std::unordered_map<node_type, uint32_t>& clusters, double etha, double mu, double epsilon)
{
      /*double etha = get_global_properties().rank_etha;
      double mu = get_global_properties().rank_mu;
      double epsilon = get_global_properties().rank_epsilon;*/

      size_t N = outlink_matrix.size1();
      size_t min_cluster = std::numeric_limits<size_t>::max();

      //condense clusters
      std::unordered_map<uint32_t, std::unordered_set<node_type>> cluster_sets;
      std::unordered_map<node_type, uint32_t> condensed_clusters;
      {
         std::unordered_map<uint32_t, std::unordered_set<node_type>> cluster_sets_sparse;
         for (const auto &p : clusters) {
            if (p.second < min_cluster) min_cluster = p.second;
            cluster_sets_sparse[p.second].insert(p.first);
         }

         size_t counter = 0;
         for (auto &p : cluster_sets_sparse) {
            auto next_cluster = counter;
            for (auto node_id : p.second)
               condensed_clusters[node_id] = next_cluster;
            cluster_sets[next_cluster] = std::move(p.second);
            ++counter;
         }
      }


      size_t clusters_size = cluster_sets.size();

      matrix_type matrixA(clusters_size, N, 0);
      for (const auto& p : condensed_clusters) {
            matrixA(p.second, p.first) = 1;
      }

      matrix_type teleportation_matrix(N, N, 0);
      double block_multiplier = 1. / clusters_size;
      for (size_t j = 0; j < N; ++j) {
            double val = block_multiplier / cluster_sets.at(condensed_clusters.at(j)).size();//cluster_sets.at(j).size();
            for (size_t i = 0; i < N; ++i)
                  teleportation_matrix(i,j) = val;
      }

      //matrix M = R * A
      matrix_type matrixR(N, clusters_size, 0);
      for (size_t i = 0; i < N; ++i) {
            std::unordered_set<uint32_t> proximal_nodes{static_cast<unsigned int>(i)};
            for (size_t i2 = 0; i2 < N; ++i2) {
                  if (outlink_matrix(i,i2) > 0)
                        proximal_nodes.insert(i2);
            }
            std::set<uint32_t> proximity_clusters;
            for (const uint32_t node : proximal_nodes) {
                  proximity_clusters.insert(condensed_clusters.at(node));
            }
            for (size_t j = 0; j < clusters_size; ++j) {
                  uint32_t cluster = j;
                  if (proximity_clusters.find(cluster) != proximity_clusters.end()) {
                        matrixR(i,j) = 1. / proximity_clusters.size() / cluster_sets.at(cluster).size();
                  }
            }
      }
      
      /*for (size_t j = 0; j < clusters_size; ++j) {
            const auto& j_set = cluster_sets.at(j);
            double val = 1 / j_set.size();
            uint32_t proximal_blocks_count = 0;
            for (size_t i = 0; i < N; ++i) {

            }
            
      }*/

      //std::cerr << "R: " << matrixR << std::endl;

      using boost::numeric::ublas::prod;
      using boost::numeric::ublas::norm_1;
      using boost::numeric::ublas::trans;

      boost::numeric::ublas::vector<double> NCDARank(N, 1), vec(N, 0); //init value!
      do{
            vec = NCDARank;
            boost::numeric::ublas::vector<double> rank_r_prod = prod(vec, matrixR);
            NCDARank = prod(etha * vec, outlink_matrix) + prod(mu * rank_r_prod, matrixA) + prod((1 - etha - mu) * vec, teleportation_matrix);
            //std::cerr << "rank: " << NCDARank << std::endl;
      } while (norm_1(NCDARank - vec) > epsilon);

      
      return NCDARank;
}

void database::calculate_importance_score(fc::time_point_sec last_maintenance)
{
      using boost::numeric::ublas::element_prod;
      using boost::numeric::ublas::norm_1;

      const auto& globs = get_global_properties().parameters;

      auto min_transfer = globs.min_transfer_for_clustering.value;
      double c_mu       = globs.clustering_mu;
      double c_epsilon  = globs.clustering_epsilon;

      double r_mu       = globs.rank_mu;
      double r_etha     = globs.rank_etha;
      double r_epsilon  = globs.rank_epsilon;

      double s_cluster  = globs.structure_cluster_weight;
      double s_outlier  = globs.structure_outlier_weight;

      double i_omega_o  = globs.importance_omega_o;
      double i_omega_i  = globs.importance_omega_i;

      std::unordered_map<uint64_t, uint64_t> accounts_map;
      auto outlink = construct_outlink_matrix(last_maintenance, accounts_map);

      boost::numeric::ublas::vector<double> outlink_vec(outlink.size1());
      for (size_t i = 0; i < outlink.size1(); ++i) {
            boost::numeric::ublas::matrix_row<outlink_matrix_type> row(outlink, i);
            outlink_vec(i) = std::accumulate(row.begin(), row.end(), 0.);
      }
      using boost::numeric::ublas::norm_1;
      for (size_t j = 0; j < outlink.size2(); ++j) {
         boost::numeric::ublas::matrix_column<outlink_matrix_type> col(outlink, j);
         double norm = norm_1(col);
         if (norm > 0)
            for (size_t i = 0; i < col.size(); ++i)
                  col(i) = col(i) / norm;
      }

      outlink = boost::numeric::ublas::trans(outlink);

      auto clusters = detail::cluster_graph_simple(outlink, accounts_map, c_epsilon, c_mu, min_transfer);
      auto rank = detail::getNCDAwareRank(outlink, clusters, r_etha, r_mu, r_epsilon);
      //std::cerr << "R " << rank << '\n';

      boost::numeric::ublas::vector<uint64_t> balances(rank.size());
      boost::numeric::ublas::vector<double> structure_weight(rank.size());
      for (const auto& p : accounts_map) {
            account_id_type acc_id(p.first);
            balances(p.second) = get_balance(acc_id, asset_id_type()).amount.value;
            
            structure_weight(p.second) = (clusters.at(p.second) > MIN_CLUSTER) ? s_cluster : s_outlier ;
      }

      /*for (size_t i = 0; i < outlink.size1(); ++i) {
            boost::numeric::ublas::matrix_column<outlink_matrix_type> col(outlink, i);
            outlink_vec(i) = std::accumulate(col.begin(), col.end(), 0.);
      }*/

      //std::cerr << "Bs " << balances << "\n";
      //std::cerr << "Ov " << outlink_vec << "\n";

      boost::numeric::ublas::vector<double> inter_vec = balances + i_omega_o * outlink_vec;
      double norm = norm_1(inter_vec);
      //std::cerr << "iv " << inter_vec / norm << "\n";
      boost::numeric::ublas::vector<double> score = element_prod(inter_vec / norm + i_omega_i * rank, structure_weight);

      //std::cerr << "Sc " << score << "\n";

      const auto& all_accounts = get_index_type<account_index>().indices().get<by_id>();
      for (auto& account : all_accounts) {
   
         //auto& sender_statistics = account.statistics(*this);
   
         //std::cout << sender_account.name << " : " << uint64_t(sender_account.id) << " = " << acc_transfers << '\n';
   
         double importance_score = 0;
         auto acc_num = object_id_type(account.get_id()).instance();
         if (accounts_map.find(acc_num) != accounts_map.end())
            importance_score = score(accounts_map.at(acc_num));

         modify(account.statistics(*this), [&](account_statistics_object& s)
         {                
            s.importance_score = importance_score;
         });
      }
}

} }