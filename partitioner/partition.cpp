#include "general_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "csrmat.hpp"
#include "numa_configuration.hpp"

#include <vector>

std::vector<std::vector<MeshMetaData>> partition(const MeshMetaData& mesh_meta,
                                                 const int num_partitions,
                                                 const int num_nodes,
                                                 const int ranks_per_locality,
                                                 const NumaConfiguration& numa_config) {
    // To do: add an additional layer of support for assigning submeshes to NUMA
    // domains
    // const int num_localities = num_nodes *
    // numa_config.get_num_numa_domains();

    std::unordered_map<int, double> element_weights;
    std::unordered_map<std::pair<int, int>, double> edge_weights;

    for (const auto& elt : mesh_meta.elements) {
        element_weights.insert(std::make_pair(elt.first, 1.));

        for (const uint neigh_id : elt.second.neighbor_ID) {
            if (neigh_id != DEFAULT_ID) {
                std::pair<int, int> edg = {std::min(elt.first, neigh_id), std::max(elt.first, neigh_id)};

                edge_weights[edg] = 1;
            }
        }
    }

    CSRMat<> mesh_graph(element_weights, edge_weights);
    std::vector<std::function<double(int)>> cons;
    cons.push_back([&element_weights](int i) { return element_weights.at(i); });

    std::vector<int64_t> mesh_part = metis_part(mesh_graph, num_partitions, cons, 1.05);

    std::unordered_map<int, int64_t> elt2partition;
    const std::vector<int>& elts_sorted = mesh_graph.node_ID();
    for (uint i = 0; i < elts_sorted.size(); ++i) {
        elt2partition.insert({elts_sorted.at(i), mesh_part.at(i)});
    }

    {
        double inter_submesh_edge_cuts(0);
        for (const auto& edg : edge_weights) {
            if (elt2partition.at(edg.first.first) != elt2partition.at(edg.first.second)) {
                inter_submesh_edge_cuts += edg.second;
            }
        }

        std::cout << "Results:\n";
        std::cout << "  Percentage of inter-submesh edge cuts: " << inter_submesh_edge_cuts / edge_weights.size() * 100
                  << " %\n";
    }

    // partition submeshes onto nodes
    std::unordered_map<int, double> submesh_weight;
    for (auto& elt : element_weights) {
        if (submesh_weight.count(elt2partition.at(elt.first))) {
            submesh_weight[elt2partition[elt.first]] += elt.second;
        } else {
            submesh_weight.insert(std::make_pair(elt2partition[elt.first], elt.second));
        }
    }

    std::unordered_map<std::pair<int, int>, double> submesh_edge_weight;
    for (auto& ew : edge_weights) {
        std::pair<int, int> sbmsh_pair{std::min(elt2partition.at(ew.first.first), elt2partition.at(ew.first.second)),
                                       std::max(elt2partition.at(ew.first.first), elt2partition.at(ew.first.second))};

        if (submesh_edge_weight.count(sbmsh_pair)) {
            submesh_edge_weight[sbmsh_pair] += ew.second;
        } else {
            submesh_edge_weight.insert(std::make_pair(sbmsh_pair, ew.second));
        }
    }

    // sbmsh_wght takes the submeshes and distribtutes them across localities
    CSRMat<> sbmsh_graph(submesh_weight, submesh_edge_weight);
    std::vector<std::function<double(int)>> sbmsh_cons;
    sbmsh_cons.push_back([&submesh_weight](int i) { return submesh_weight.at(i); });

    std::vector<int64_t> sbmsh_part = metis_part(sbmsh_graph, num_nodes, sbmsh_cons, 1.02);

    std::vector<int> permutation;
    permutation.reserve(sbmsh_part.size());
    for (uint i = 0; i < sbmsh_part.size(); ++i) {
        permutation.push_back(i);
    }

    std::unordered_map<int, int64_t> partition2node;
    std::unordered_map<int, int64_t> partition2rank;
    const std::vector<int>& sbmshs_sorted = sbmsh_graph.node_ID();
    for (uint i = 0; i < sbmshs_sorted.size(); ++i) {
        partition2node.insert({sbmshs_sorted.at(i), sbmsh_part.at(i)});
    }

    std::vector<std::vector<MeshMetaData>> submeshes(num_nodes * ranks_per_locality);
    std::vector<uint> partition2local_partition(num_partitions);
    {
        std::vector<uint> local_partition_counter(num_nodes * ranks_per_locality, 0);
        std::vector<uint> local_rank_counter(num_nodes, 0);
        for (auto& p_n : partition2node) {
            // Assign submeshes to ranks on a node, by cycling round robin through the meshes
            const uint rank = local_rank_counter[p_n.second] + ranks_per_locality * p_n.second;
            local_rank_counter[p_n.second] = (local_rank_counter[p_n.second] + 1) % ranks_per_locality;

            partition2rank[p_n.first] = rank;
            partition2local_partition[p_n.first] = local_partition_counter[rank]++;
        }

        for (uint sm = 0; sm < submeshes.size(); ++sm) {
            submeshes[sm].resize(local_partition_counter[sm]);
            for (uint i = 0; i < local_partition_counter[sm]; ++i) {
                submeshes[sm][i].mesh_name = mesh_meta.mesh_name + "_" + std::to_string(sm) + "_" + std::to_string(i);
            }
        }
    }

    {  // assemble submeshes
        for (const auto& elt : mesh_meta.elements) {
            uint partition = elt2partition.at(elt.first);
            uint rank = partition2rank[partition];
            uint loc_part = partition2local_partition[partition];

            // error here the partition in question is not actually the
            // partition you are interested in.
            submeshes[rank][loc_part].elements.insert(elt);

            for (uint id : elt.second.node_ID) {
                submeshes[rank][loc_part].nodes[id] = mesh_meta.nodes.at(id);
            }
        }
    }

    {  // update boundaries for split edges
        for (const auto& edg : edge_weights) {
            if (elt2partition.at(edg.first.first) != elt2partition.at(edg.first.second)) {
                uint elt_A = edg.first.first;
                uint part_A = elt2partition.at(elt_A);
                uint loc_part_A = partition2local_partition[part_A];
                uint rank_A = partition2rank[part_A];

                uint elt_B = edg.first.second;
                uint part_B = elt2partition.at(elt_B);
                uint loc_part_B = partition2local_partition[part_B];
                uint rank_B = partition2rank[part_B];

                for (uint k = 0; k < 3; ++k) {
                    if (elt_B == mesh_meta.elements.at(elt_A).neighbor_ID[k]) {
                        ElementMetaData& curr_elt = submeshes[rank_A][loc_part_A].elements.at(elt_A);
                        curr_elt.neighbor_ID[k] = DEFAULT_ID;
                        curr_elt.boundary_type[k] = SWE::BoundaryConditions::distributed;
                    }
                }

                for (uint k = 0; k < 3; ++k) {
                    if (elt_A == mesh_meta.elements.at(elt_B).neighbor_ID[k]) {
                        ElementMetaData& curr_elt = submeshes[rank_B][loc_part_B].elements.at(elt_B);
                        curr_elt.neighbor_ID[k] = DEFAULT_ID;
                        curr_elt.boundary_type[k] = SWE::BoundaryConditions::distributed;
                    }
                }
            }
        }
    }

    {  // compute percentage of inter-node edge cuts
        double inter_node_edge_cuts(0);
        for (const auto& edg : submesh_edge_weight) {
            if (partition2node.at(edg.first.first) != partition2node.at(edg.first.second)) {
                inter_node_edge_cuts += edg.second;
            }
        }

        std::cout << "  Percentage of inter-node edge cuts: " << inter_node_edge_cuts / edge_weights.size() * 100
                  << " %\n";
    }

    {  // compute imbalance across nodes
        double max_load(0);
        double avg_load(0);
        std::vector<double> node_weight(num_nodes, 0);

        for (const auto& sw : submesh_weight) {
            node_weight.at(partition2node.at(sw.first)) += sw.second;
        }

        for (const auto& wght : node_weight) {
            if (wght > max_load) {
                max_load = wght;
            }
            avg_load += wght;
        }

        avg_load /= num_nodes;
        double imbalance = (max_load - avg_load) / avg_load;

        std::cout << "  Imbalance across NUMA domains: " << imbalance << '\n';
    }

    return submeshes;
}
