#include "general_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include "csrmat.hpp"
#include "partition_type.hpp"

#include <memory>
#include <vector>

void refine_strict_tracked(const CSRMat<>& g,
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
    assert( submesh_graph.size() == num_partitions );

    std::vector<std::function<double(int)>> cons2;
    cons2.push_back([&submesh_graph](int ID){ return submesh_graph.node_weight(ID); });

    uint num_ranks = ranks_per_numa_domain * numa_domains_per_node * num_nodes;
    std::vector<int64_t> rank_partition_vec = metis_part(submesh_graph, num_ranks, cons2, 1.02);

    std::vector<PartitionType> submesh_partitions;
    submesh_partitions.emplace_back(num_ranks,
                                    submesh_graph,
                                    rank_partition_vec);

    for ( auto& sp : submesh_partitions ) {
      std::cout << sp << std::endl;
    }

    uint submeshes_per_rank = num_partitions / num_ranks;
    std::cout << "submeshes_per_rank = " << submeshes_per_rank << '\n';
    assert( submeshes_per_rank * num_ranks == num_partitions );

    std::vector<uint> coarsening_factors{submeshes_per_rank * ranks_per_numa_domain, ranks_per_numa_domain, numa_domains_per_node};
    std::cout << "  coarsening_factors = { " << coarsening_factors[0] << ", " << coarsening_factors[1] << ", " << coarsening_factors[2] << "}" << std::endl;
    std::cout << "submesh_partitions.size() = " << submesh_partitions.size() << '\n';

    refine_strict_tracked(submesh_graph, submesh_partitions, coarsening_factors);

    std::vector<std::vector<MeshMetaData>> submeshes(num_ranks);

    std::unordered_map<int64_t,int64_t> rank2local_rank;
    std::unordered_map<int64_t,int64_t> numa2local_numa;
    std::unordered_map<int64_t,int64_t> part2local_part;
    {
      std::vector<int64_t> local_part_counter(num_ranks,0);
        std::vector<int64_t> local_rank_counter(numa_domains_per_node*num_nodes,0);
        std::vector<int64_t> local_numa_counter(num_nodes,0);

        for ( auto& v_p : submesh_partitions[RANK].vertex2partition ) {
          int64_t part_id = v_p.first;
          int64_t rank_id = v_p.second;

          if ( !rank2local_rank.count(rank_id) ) {

            int64_t node_id = submesh_partitions[NODE].vertex2partition.at(part_id);
            int64_t numa_id = submesh_partitions[NUMA].vertex2partition.at(part_id);

            if ( !numa2local_numa.count(numa_id) ) {
              numa2local_numa[numa_id] = local_numa_counter[node_id]++;
            }

            int64_t flat_numa_id = node_id*numa_domains_per_node + numa2local_numa[numa_id];

            rank2local_rank[rank_id] = local_rank_counter[flat_numa_id]++;

          }
          part2local_part[part_id] = local_part_counter[rank_id]++;
        }

        for (uint sm = 0; sm < submeshes.size(); ++sm) {
            submeshes[sm].resize(local_part_counter[sm]);
            for (uint i = 0; i < local_part_counter[sm]; ++i) {
                submeshes[sm][i].mesh_name = mesh_meta.mesh_name + "_" + std::to_string(sm) + "_" + std::to_string(i);
            }
        }
    }

    {  // assemble submeshes
        for (const auto& elt : mesh_meta.elements) {
            uint partition = mesh_partition.vertex2partition.at(elt.first);
            uint numa_id = submesh_partitions[NUMA].vertex2partition.at(partition);
            uint rank_id = submesh_partitions[RANK].vertex2partition.at(partition);

            uint rank = (submesh_partitions[NODE].vertex2partition.at(partition)*numa_domains_per_node + numa2local_numa.at(numa_id))
              *ranks_per_numa_domain + rank2local_rank.at(rank_id);

            uint loc_part = part2local_part.at(partition);

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
                uint rank_id_A = submesh_partitions[RANK].vertex2partition.at(partition_A);
                uint rank_A = (numa_domains_per_node * submesh_partitions[NODE].vertex2partition.at(partition_A)
                               + numa2local_numa.at(numa_id_A))*ranks_per_numa_domain + rank2local_rank.at(rank_id_A);
                uint loc_part_A = part2local_part[partition_A];

                uint elt_B = edg.first.second;
                uint partition_B = mesh_partition.vertex2partition.at(elt_B);
                uint numa_id_B = submesh_partitions[NUMA].vertex2partition.at(partition_B);
                uint rank_id_B = submesh_partitions[RANK].vertex2partition.at(partition_B);
                uint rank_B = (numa_domains_per_node * submesh_partitions[NODE].vertex2partition.at(partition_B)
                               + numa2local_numa[numa_id_B])*ranks_per_numa_domain + rank2local_rank.at(rank_id_B);;
                uint loc_part_B = part2local_part[partition_B];


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