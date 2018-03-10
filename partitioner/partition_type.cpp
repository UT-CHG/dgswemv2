#include "partition_type.hpp"

PartitionType::PartitionType(int64_t num_partitions, const CSRMat<>& g, const std::vector<int64_t>& partition) : num_partitions(num_partitions) {
    for ( uint indx = 0; indx < g.size(); ++indx ) {
        vertex2partition.insert(std::make_pair(g.node_ids().at(indx), partition[indx]));
    }
}

CSRMat<> PartitionType::make_partition_graph(const CSRMat<>& g) {
    std::unordered_map<int,double> super_nodes;
    std::unordered_map<std::pair<int,int>, double> super_edges;

    for ( auto& n_w : g.get_node_wghts_map() ) {
        if ( vertex2partition.count(n_w.first) ) {
            int partition = vertex2partition[n_w.first];
            if ( super_nodes.count(partition) ) {
                super_nodes[partition] += n_w.second;
            } else {
                super_nodes[partition] = n_w.second;
            }
        }
    }

    for ( auto& e_w : g.get_edge_wgts_map() ) {
        if ( vertex2partition.count(e_w.first.first) &&
             vertex2partition.count(e_w.first.second) ) {
            std::pair<int,int> super_edge_name{ std::min( vertex2partition.at(e_w.first.first),
                                                          vertex2partition.at(e_w.first.second) ),
                                                std::max( vertex2partition.at(e_w.first.first),
                                                          vertex2partition.at(e_w.first.second) ) };
            if ( super_edges.count( super_edge_name ) ) {
                super_edges[super_edge_name] += e_w.second;
            } else {
                super_edges[super_edge_name] = e_w.second;
            }
        }
    }

    return CSRMat<>(super_nodes, super_edges);
}