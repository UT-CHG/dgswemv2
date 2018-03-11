#include "partition_type.hpp"
#include "csrmat.hpp"
#include "steal_queue.hpp"

#include <deque>

void manually_refine2(const CSRMat<>& g, PartitionType&p, bool override=false) {
//by convention here, gain is the change in edge cuts, but going from 0 to 1
    std::array<uint,2> discrete_weights{0,0};
    std::array<std::deque<int>,2> elements;

    for ( auto& v_p : p.vertex2partition ) {
        discrete_weights[ v_p.second ] += 1;
        elements[ v_p.second ].push_back(v_p.first);
    }

    if ( discrete_weights[0] == discrete_weights[1] ) {
        return;
    }

    if ( ( (discrete_weights[0] + discrete_weights[1]) % 2 != 0) && !override ) {
        std::string err_msg{"Error: Unable to exactly manually_refine2 partitions with weights: "
                + std::to_string(discrete_weights[0])+", "+std::to_string(discrete_weights[1])+"\n"};
        throw std::logic_error(err_msg);
    }

    uint high = discrete_weights[1] > discrete_weights[0];
    uint low = 1 - high;

    std::cout << "high: " << high << " low: " << low << '\n';

    //In the case of an override there will be an odd number of tiles, and we would like
    // to ensure that one extra tile ends up on the low end. Note this number is strictly positive.
    uint num_tiles_to_be_moved = (discrete_weights[high] - discrete_weights[low] + 1)/2;
    std::cout << "High elements:";
    for( auto& elt : elements[high] ) {
        std::cout << ' ' << elt;
    }
    std::cout << '\n';

    StealQ steal_q(g, p, elements[high]);

    std::cout << "num_tiles_to_be_moved: " << num_tiles_to_be_moved << '\n';

    while ( num_tiles_to_be_moved > 0 ) {
        int id = steal_q.steal_one();
        std::cout << "stealing id: " << id << '\n';
        assert( p.vertex2partition[id] == high );
        p.vertex2partition[id] = low;
        --num_tiles_to_be_moved;
    }
}
/*
void manually_refine3() {

}


*/
PartitionType coarsen(const CSRMat<>& g, const PartitionType& p, uint coarsening_factor) {
    //if exact these are the number of partitions
    int64_t num_partitions = p.vertex2partition.size()/coarsening_factor;

    //Otherwise, add an extra partition
    if ( num_partitions*coarsening_factor != p.num_partitions ) {
        num_partitions = static_cast<int64_t>( std::ceil( static_cast<double>(g.size())/coarsening_factor));
    }

    CSRMat<> g_super = p.make_partition_graph();

    std::vector<std::function<double(int)>> cons;
    cons.push_back([&g_super](int id) { return g_super.node_weight(id); });

    std::vector<int64_t> coarse_partition = metis_part(g_super, num_partitions, cons, 1.02);

    return PartitionType(num_partitions,
                         g, p,
                         coarse_partition);
}

void refine_strict(const CSRMat<>& g, PartitionType& p, uint coarsening_factor) {
    if ( coarsening_factor <= 1 ) {
        throw std::logic_error("Error in refine_strict;"
                               "Coarsening factor needs to be greater than 1\n");
    }

    //Take care of base cases
    if ( p.num_partitions == 1 ) {
        return;
    } else if ( p.num_partitions == 2 ) {
        manually_refine2(g,p);
        return;
/*    } else if ( p.num_partitions ==3 ) {
        manually_refine3(g,p);
        return;*/
    }

    //Nothing to be done
    //if ( p.is_balanced(g,p) ) {
    //    return;
    //}

    PartitionType p_coarse = coarsen(g, p, coarsening_factor);

/*
    //save copy of old vertex placements to track movement
    PartitionType p_coarse_old = p_coarse;
    refine_strict(g, p_coarse, 2);

    std::unordered_set<int> moved_vertices;
    for ( auto& v_p : p_coarse.vertex2partition ) {
        if ( v_p.second != p_coarse_old.vertex2partition.at(v_p.first) ) {
            moved_vertices.insert(v_p.first);
            //moved vertices are temporarily not assinged to any partition
            p[v_p.first] = -1;
        }
    }

    //make sub partitioned subgraphs
    std::vector<PartitionType> sub_partitions(p_coarse.num_partitions);
    {
        std::vector<std::unordered_set<int64_t>> partitions_on_coarse_partition(p_coarse.num_partitions);
        for ( auto& v_p : p_coarse.vertex2partition ) {
            if ( p[v_p.first] != -1 ) {
                int64_t fine_partition = p[v_p.first];
                sub_partitions.vertex2partition[v_p.first] = fine_partition;
                partitions_on_coarse_partition[v_p.second].insert(fine_partition);
            }
        }

        }*/

    
}