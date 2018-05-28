#include "general_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "csrmat.hpp"

#include <vector>

std::vector<std::vector<MeshMetaData>> partition(const MeshMetaData& mesh_meta,
                                                 const std::unordered_map<int, std::vector<double>>& problem_weights,
                                                 const int num_partitions,
                                                 const int num_nodes,
                                                 const int ranks_per_locality,
                                                 const bool rank_balanced) {
    std::unordered_map<int, std::vector<double>> element_weights;
    std::unordered_map<std::pair<int, int>, double> edge_weights;

    for (const auto& elt : mesh_meta.elements) {
        if (!rank_balanced) {
            element_weights.insert(std::make_pair(elt.first, std::vector<double>({1.})));
        }

        for (const uint neigh_id : elt.second.neighbor_ID) {
            if (neigh_id != DEFAULT_ID) {
                std::pair<int, int> edg = {std::min(elt.first, neigh_id), std::max(elt.first, neigh_id)};

                edge_weights[edg] = 1;
            }
        }
    }

    CSRMat mesh_graph(rank_balanced ? problem_weights : element_weights, edge_weights);

    std::vector<int64_t> mesh_part = metis_part(mesh_graph, num_partitions, 1.05);

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
    std::unordered_map<int, std::vector<double>> submesh_weight;
    for (auto& elt : problem_weights) {
        if (submesh_weight.count(elt2partition.at(elt.first))) {
            for (uint c = 0; c < elt.second.size(); ++c) {
                submesh_weight[elt2partition[elt.first]][c] += elt.second[c];
            }
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
    CSRMat sbmsh_graph(submesh_weight, submesh_edge_weight);

    std::vector<int64_t> sbmsh_part = metis_part(sbmsh_graph, num_nodes, 1.02);

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
            const uint rank                = local_rank_counter[p_n.second] + ranks_per_locality * p_n.second;
            local_rank_counter[p_n.second] = (local_rank_counter[p_n.second] + 1) % ranks_per_locality;

            partition2rank[p_n.first]            = rank;
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
            uint rank      = partition2rank[partition];
            uint loc_part  = partition2local_partition[partition];

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
                uint elt_A      = edg.first.first;
                uint part_A     = elt2partition.at(elt_A);
                uint loc_part_A = partition2local_partition[part_A];
                uint rank_A     = partition2rank[part_A];

                uint elt_B      = edg.first.second;
                uint part_B     = elt2partition.at(elt_B);
                uint loc_part_B = partition2local_partition[part_B];
                uint rank_B     = partition2rank[part_B];

                for (uint k = 0; k < 3; ++k) {
                    if (elt_B == mesh_meta.elements.at(elt_A).neighbor_ID[k]) {
                        ElementMetaData& curr_elt = submeshes[rank_A][loc_part_A].elements.at(elt_A);
                        curr_elt.neighbor_ID[k]   = DEFAULT_ID;
                        curr_elt.boundary_type[k] = SWE::BoundaryTypes::distributed;
                    }
                }

                for (uint k = 0; k < 3; ++k) {
                    if (elt_A == mesh_meta.elements.at(elt_B).neighbor_ID[k]) {
                        ElementMetaData& curr_elt = submeshes[rank_B][loc_part_B].elements.at(elt_B);
                        curr_elt.neighbor_ID[k]   = DEFAULT_ID;
                        curr_elt.boundary_type[k] = SWE::BoundaryTypes::distributed;
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
        std::size_t num_constraints = submesh_weight.cbegin()->second.size();
        std::vector<double> max_load(num_constraints, 0);
        std::vector<double> avg_load(num_constraints, 0);

        std::vector<std::vector<double>> node_weight(num_nodes, std::vector<double>(num_constraints, 0));

        for (const auto& sw : submesh_weight) {
            for (uint c = 0; c < num_constraints; ++c) {
                node_weight.at(partition2node.at(sw.first))[c] += sw.second[c];
            }
        }

        for (const auto& wght : node_weight) {
            for (uint c = 0; c < num_constraints; ++c) {
                max_load[c] = std::max(wght[c], max_load[c]);
                avg_load[c] += wght[c];
            }
        }

        std::cout << "  Imbalance across NUMA domains:\n";
        for (uint c = 0; c < num_constraints; ++c) {
            avg_load[c] /= num_nodes;
            double imbalance = (max_load[c] - avg_load[c]) / avg_load[c];
            std::cout << "    Constraint " << c << ": " << imbalance << '\n';
        }
        std::cout << '\n';
    }

    return submeshes;
}
