#include "general_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "csrmat.hpp"
#include "partition_type.hpp"

#include <vector>

std::vector<PartitionType> refine_strict_tracked(const CSRMat<>& g,
                                                 std::vector<PartitionType>& ps,
                                                 const std::vector<uint>& c);

namespace {
constexpr uint RANK{0u};
constexpr uint NUMA{1u};
constexpr uint NODE{2u};
}

std::vector<std::vector<MeshMetaData>> partition2(const MeshMetaData& mesh_meta,
                                                 const int num_partitions,
                                                 const uint ranks_per_numa_domain,
                                                 const uint numa_domains_per_node,
                                                 const uint num_nodes) {

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

    //Partition mesh into submeshes
    CSRMat<> mesh_graph(element_weights, edge_weights);

    std::vector<std::function<double(int)>> cons;
    cons.push_back([&element_weights](int ID){ return element_weights.at(ID); });

    std::vector<int64_t> submesh_partition_vec = metis_part(mesh_graph, num_partitions, cons, 1.05);

    PartitionType mesh_partition(num_partitions, mesh_graph, submesh_partition_vec);

    //Partition Submeshes onto ranks
    CSRMat<> submesh_graph = mesh_partition.make_partition_graph();

    std::vector<PartitionType> submesh_partitions;
    submesh_partitions.emplace_back(ranks_per_numa_domain * numa_domains_per_node * num_nodes,
                                    mesh_graph,
                                    submesh_partition_vec);

    std::vector<uint> coarsening_factors{ranks_per_numa_domain, numa_domains_per_node, num_nodes};

    refine_strict_tracked(mesh_graph, submesh_partitions, coarsening_factors);

    std::vector<std::vector<MeshMetaData>> submeshes(ranks_per_numa_domain * numa_domains_per_node * num_nodes);

    std::unordered_map<int64_t,int64_t> partition2local_partition;
    std::unordered_map<int64_t,int64_t> numa2local_numa;
    {
        std::vector<int64_t> local_partition_counter(numa_domains_per_node*num_nodes,0);
        std::vector<int64_t> local_numa_counter(num_nodes,0);

        for ( auto& v_p : submesh_partitions[RANK].vertex2partition ) {
            int64_t rank_id = v_p.second;
            if ( !partition2local_partition.count(rank_id) ) {
                int submesh_id = v_p.first;

                int64_t node_id = submesh_partitions[NODE].vertex2partition.at(submesh_id);
                int64_t numa_id = submesh_partitions[NUMA].vertex2partition.at(submesh_id);

                if ( !numa2local_numa.count(numa_id) ) {
                    numa2local_numa[numa_id] = local_numa_counter[node_id]++;
                }

                int64_t flat_numa_id = node_id*numa_domains_per_node + numa2local_numa[numa_id];

                partition2local_partition[rank_id] = local_partition_counter[flat_numa_id]++;

            }
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
            uint partition = mesh_partition.vertex2partition.at(elt.first);
            uint numa_id = submesh_partitions[NUMA].vertex2partition.at(partition);
            uint rank = numa_domains_per_node * submesh_partitions[NODE].vertex2partition.at(partition)
                + numa2local_numa[numa_id];

            uint loc_part = partition2local_partition[partition];

            // error here the partition in question is not actually the
            // partition you are interested in.
            submeshes[rank][loc_part].elements.insert(elt);

            for (uint id : elt.second.node_ids) {
                submeshes[rank][loc_part].nodes[id] = mesh_meta.nodes.at(id);
            }
        }
    }

    {  // update boundaries for split edges
        for (const auto& edg : edge_weights) {
            if (mesh_partition.vertex2partition.at(edg.first.first) != mesh_partition.vertex2partition.at(edg.first.second)) {
                uint elt_A = edg.first.first;
                uint partition_A = mesh_partition.vertex2partition.at(elt_A);
                uint numa_id_A = submesh_partitions[NUMA].vertex2partition.at(partition_A);
                uint rank_A = numa_domains_per_node * submesh_partitions[NODE].vertex2partition.at(partition_A)
                    + numa2local_numa[numa_id_A];
                uint loc_part_A = partition2local_partition[partition_A];

                uint elt_B = edg.first.second;
                uint partition_B = mesh_partition.vertex2partition.at(elt_B);
                uint numa_id_B = submesh_partitions[NUMA].vertex2partition.at(partition_B);
                uint rank_B = numa_domains_per_node * submesh_partitions[NODE].vertex2partition.at(partition_B)
                    + numa2local_numa[numa_id_B];
                uint loc_part_B = partition2local_partition[partition_B];


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

    return submeshes;
}