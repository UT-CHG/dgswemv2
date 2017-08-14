#include "general_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "csrmat.hpp"
#include "numa_configuration.hpp"

#include <vector>

std::vector<std::vector<MeshMetaData> > partition(const MeshMetaData& mesh_meta,
                                                  const int num_partitions,
                                                  const int num_nodes,
                                                  const NumaConfiguration& numa_config)
{
    const int num_localities = num_nodes*numa_config.get_num_numa_domains();

    std::unordered_map<int, double> element_weights;
    std::unordered_map<std::pair<int,int>, double> edge_weights;

    for ( const auto& elt : mesh_meta._elements ) {
        element_weights.insert( std::make_pair( elt.first, 1.) );

        for ( const uint neigh_id : elt.second.neighbor_ID ) {
            if ( neigh_id != DEFAULT_ID ) {
                std::pair<int,int> edg = { std::min( elt.first, neigh_id),
                                           std::max( elt.first, neigh_id) };

                edge_weights[edg] = 1;
            }
        }
    }


    CSRMat<> mesh_graph(element_weights, edge_weights);
    std::vector<std::function<double(int)> > cons;
    cons.push_back([&element_weights](int i){
            return element_weights.at(i);});

    std::vector<int64_t> mesh_part = metis_part( mesh_graph, num_partitions,
                                                 cons, 1.05 );

    std::unordered_map<int,int64_t> elt2partition;
    const std::vector<int>& elts_sorted = mesh_graph.node_ids();
    for ( uint i = 0; i < elts_sorted.size(); ++i ) {
        elt2partition.insert( { elts_sorted.at(i), mesh_part.at(i) } );
    }

    {
        double inter_submesh_edge_cuts(0);
        for ( const auto& edg : edge_weights ) {
            if ( elt2partition.at(edg.first.first) != elt2partition.at(edg.first.second) ) {
                inter_submesh_edge_cuts += edg.second;
            }
        }

        std::cout << "  Percentage of inter-submesh edge cuts: "
                  << inter_submesh_edge_cuts/edge_weights.size()*100
                  << " %\n";
    }

    /*//partition submeshes onto localities
    std::unordered_map<int,double> submesh_weight;
    for ( auto& elt : element_weights ) {
        if ( submesh_weight.count( elt2partition.at(elt.first) ) ) {
            submesh_weight[ elt2partition[ elt.first ] ] += elt.second;
        } else {
            submesh_weight.insert( std::make_pair( elt2partition[elt.first], elt.second ));
        }
    }

    std::unordered_map<std::pair<int,int>,double> submesh_edge_weight;
    for (auto& ew : edge_weights ) {
        std::pair<int,int> sbmsh_pair { std::min( elt2partition.at(ew.first.first ),
                                                      elt2partition.at(ew.first.second)),
                std::max( elt2partition.at(ew.first.first ),
                          elt2partition.at(ew.first.second)) };

        if ( submesh_edge_weight.count( sbmsh_pair ) ) {
            submesh_edge_weight[sbmsh_pair] += ew.second;
        } else {
            submesh_edge_weight.insert( std::make_pair( sbmsh_pair, ew.second ) );
        }
    }

    //sbmsh_wght takes the submeshes and distribtutes them across localities
    CSRMat<> sbmsh_graph(submesh_weight,submesh_edge_weight);
    std::vector<std::function<double(int)> > sbmsh_cons;
    sbmsh_cons.push_back([&submesh_weight](int i){
            return submesh_weight.at(i);});


    std::vector<int64_t> sbmsh_part = metis_part( sbmsh_graph, num_localities,
                                                  sbmsh_cons, 1.02 );

    std::vector<int> permutation;
    permutation.reserve(sbmsh_part.size());
    for ( uint i = 0; i < sbmsh_part.size(); ++i ) {
        permutation.push_back(i);
    }

    std::unordered_map<int,int64_t> partition2node;
    const std::vector<int>& sbmshs_sorted = sbmsh_graph.node_ids();
    for ( uint i = 0; i < sbmshs_sorted.size(); ++i ) {
        partition2node.insert( { sbmshs_sorted.at(i), sbmsh_part.at(i) } );
    }


    { // compute percentage of inter-node edge cuts
        double inter_node_edge_cuts(0);
        for ( const auto& edg : submesh_edge_weight ) {
            if ( partition2node.at(edg.first.first) != partition2node.at(edg.first.second) ) {
                inter_node_edge_cuts += edg.second;
            }
        }

        std::cout << "  Percentage of inter-node edge cuts: "
                  << inter_node_edge_cuts/edge_weights.size()*100
                  << " %\n";

    }

    { //compute imbalance across NUMA Domains
        double max_load(0);
        double avg_load(0);
        for ( const auto& sw : submesh_weight ) {
            if ( sw.second > max_load ) {
                max_load = sw.second;
            }
            avg_load += sw.second;
        }
        avg_load /= submesh_weight.size();
        double imbalance = ( max_load - avg_load ) / avg_load;

        std::cout << "  Imbalance across NUMA domains: " << imbalance << '\n';
        }*/

    return std::vector<std::vector<MeshMetaData>>();
}
