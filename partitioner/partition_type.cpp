#include "partition_type.hpp"
#include "utilities/almost_equal.hpp"

PartitionType::PartitionType(int64_t num_partitions,
                             const CSRMat<>& g,
                             const std::vector<int64_t>& partition)
    : num_partitions(num_partitions), g_ref(g), ideal_weights(num_partitions,1)
{
    assert( num_partitions > 0 );
    assert( g.size() == partition.size() );

    for ( uint id = 0; id < g.size(); ++id ) {
        vertex2partition.insert( {g.node_ids()[id], partition[id]} );
    }
}

PartitionType::PartitionType(int64_t num_partitions,
                             const CSRMat<>& g,
                             const std::vector<int64_t>& partition,
                             const std::vector<double>& ideal_weights_)
    : PartitionType(num_partitions,g,partition)
{
    ideal_weights = ideal_weights_;
}

PartitionType::PartitionType(int64_t num_partitions,
                             const CSRMat<>& g,
                             const std::vector<int>& node_ids,
                             const std::vector<int64_t>& partition,
                             const std::vector<double>& ideal_weights)
    : num_partitions(num_partitions), g_ref(g), ideal_weights(ideal_weights) {

    assert( num_partitions > 0 );
    assert( node_ids.size() == partition.size() );
    assert( ideal_weights.size() == (uint)num_partitions );

    for ( uint indx = 0; indx < node_ids.size(); ++indx ) {
        vertex2partition.insert(std::make_pair(node_ids[indx], partition[indx]));
    }
}

PartitionType::PartitionType(int64_t num_partitions,
                             const CSRMat<>& g,
                             const PartitionType& p,
                             const std::vector<int>& coarse_nodes,
                             const std::vector<int64_t>& coarse_partition)
    : num_partitions(num_partitions), g_ref(g), vertex2partition(p.vertex2partition), ideal_weights(num_partitions,0) {
    std::unordered_map<int,int64_t> coarse_partition_id;
    for ( uint i =0; i < coarse_nodes.size(); ++i ) {
      coarse_partition_id.insert( {coarse_nodes[i], coarse_partition[i]} );
    }

    for ( auto& v_p : vertex2partition ) {
        v_p.second = coarse_partition_id.at(v_p.second);
    }

    for ( uint fine_partition = 0; fine_partition < p.num_partitions; ++fine_partition ) {
      ideal_weights[ coarse_partition_id.at(fine_partition) ] += p.ideal_weights.at(fine_partition);
    }
}

CSRMat<> PartitionType::make_partition_graph() const {
    std::unordered_map<int,double> super_nodes;
    std::unordered_map<std::pair<int,int>, double> super_edges;

    for ( auto& n_w : g_ref.get_node_wghts_map() ) {
        if ( vertex2partition.count(n_w.first) ) {
            int partition = vertex2partition.at(n_w.first);
            if ( super_nodes.count(partition) ) {
                super_nodes[partition] += n_w.second;
            } else {
                super_nodes[partition] = n_w.second;
            }
        }
    }

    for ( auto& e_w : g_ref.get_edge_wgts_map() ) {
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

bool PartitionType::is_balanced() {
    //easy case
    if ( num_partitions == 1 ) { return true; }

    std::vector<double> partition_weights(num_partitions, 0);

    for ( auto& v_p : vertex2partition ) {
        partition_weights[v_p.second] += 1;
    }

    for ( uint part_id = 0; part_id < num_partitions; ++part_id ) {
        partition_weights[part_id] /= ideal_weights[part_id];
    }

    for ( uint part_id = 1; part_id < num_partitions; ++part_id ) {
        if ( !Utilities::almost_equal( partition_weights[part_id], partition_weights[0]) ) {
            return false;
        }
    }

    return true;
}

void PartitionType::add_vertex(int ID) {
    assert( !vertex2partition.count(ID) );
    int64_t partition{0};

    if ( num_partitions > 1 ) {
        std::vector<double> gain(num_partitions,0);

        for ( int neigh : g_ref.get_neighbors(ID) ) {
            if ( vertex2partition.count(neigh) ) {
                std::pair<int,int> edge_name{ std::min(ID, neigh), std::max(ID,neigh) };
                gain[vertex2partition.at(neigh)] += g_ref.edge_weight(edge_name);
            }
        }

        double max_gain=gain[0];
        for ( uint i = 1; i < partition; ++i ) {
            if ( gain[i] > max_gain ) {
                max_gain = gain[i];
                partition = i;
            }
        }
    }

    vertex2partition.insert({ID,partition});
}

std::ostream& operator<<(std::ostream& s, const PartitionType& p) {
    s << "num_partitions: " << p.num_partitions << '\n'
      << "vertex2partition: " << '\n';
    for ( auto& v_p : p.vertex2partition ) {
        s << "  " << v_p.first << " : " << v_p.second << '\n';
    }

    s << "ideal_weights: " << '\n';
    for ( uint i = 0; i < p.num_partitions; ++i ) {
        s << "  " << i << " : " << p.ideal_weights[i] << '\n';
    }
    return s;
}