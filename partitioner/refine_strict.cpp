#include "partition_type.hpp"
#include "csrmat.hpp"
#include "steal_queue.hpp"

#include "utilities/almost_equal.hpp"

#include <deque>

void manually_refine2(const CSRMat<>& g, PartitionType&p,
                      const std::vector<double>& custom_weights = std::vector<double>()) {
    //custom weights should either be empty or have size 2
    assert( custom_weights.empty() || custom_weights.size() == 2 );

    const std::vector<double>& ideal_weights = custom_weights.empty() ? p.ideal_weights : custom_weights;

//by convention here, gain is the change in edge cuts, but going from 0 to 1
    std::array<uint,2> weights{0,0};
    std::array<std::deque<int>,2> elements;

    for ( auto& v_p : p.vertex2partition ) {
        weights[ v_p.second ] += 1; //enforce unit weight
        elements[ v_p.second ].push_back(v_p.first);
    }

    uint num_tiles_to_be_moved;
    { // check feasibility, by ensuring than an integral number of tiles needs to be moved
      //  to balance the 2 partitions
        double num_tiles_fp = ( ideal_weights[1]*weights[0] - ideal_weights[0]*weights[1] )/
            (ideal_weights[0] + ideal_weights[1]);
        if ( !Utilities::almost_equal(num_tiles_fp, std::round(num_tiles_fp)) ) {
            std::string err_msg{"Error: Unable to exactly manually_refine2 partitions with weights: "
                    + std::to_string(weights[0])+", "+std::to_string(weights[1])+"\n"};
            throw std::logic_error(err_msg);
        }

        num_tiles_to_be_moved = std::abs( num_tiles_fp );
    }

    uint high = weights[1]/p.ideal_weights[1] > weights[0]/p.ideal_weights[0];
    uint low = 1 - high;

    std::cout << "high: " << high << " low: " << low << '\n';

    //In the case of an override there will be an odd number of tiles, and we would like
    // to ensure that one extra tile ends up on the low end. Note this number is strictly positive.
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

    //Nothing to be done
    if ( p.is_balanced() ) {
        return;
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


    PartitionType p_coarse = coarsen(g, p, coarsening_factor);
    std::unordered_set<int> moved_vertices;


    { //save copy of old vertex placements to track movement
        PartitionType p_coarse_old = p_coarse;
        refine_strict(g, p_coarse, 2);

        for ( auto& v_p : p_coarse.vertex2partition ) {
            if ( v_p.second != p_coarse_old.vertex2partition.at(v_p.first) ) {
                moved_vertices.insert(v_p.first);
                //moved vertices are temporarily not assinged to any partition
                p.vertex2partition[v_p.first] = -1;
            }
        }
    }

    { //make sub partitioned subgraphs
        std::vector<std::pair<int64_t,int64_t>> glo2loc(p.num_partitions,std::make_pair(-1,-1));
        std::unordered_map<std::pair<int64_t,int64_t>, int64_t> loc2glo;
        std::vector<std::size_t> local_partition_counter(p_coarse.num_partitions,0);

        for ( auto& v_p : p.vertex2partition ) {
            if ( v_p.second != -1 ) { //hasn't been moved
                if ( glo2loc[v_p.second].first == -1 ) { //hasn't been assigned yet
                    int64_t coarse_partition = p_coarse.vertex2partition[v_p.first];
                    int64_t local_partition = local_partition_counter[coarse_partition]++;

                    glo2loc[v_p.second].first = coarse_partition;
                    glo2loc[v_p.second].second = local_partition;

                    std::pair<int64_t,int64_t> local_name{coarse_partition, local_partition};
                    loc2glo[local_name] = v_p.second;
                }
            }
        }

        //build the local sub_partitions
        std::vector<PartitionType> sub_partitions;
        {
            sub_partitions.reserve(p_coarse.num_partitions);

            std::vector<std::vector<int>> node_ids(p_coarse.num_partitions);
            std::vector<std::vector<int64_t>> partition(p_coarse.num_partitions);
            for ( auto& v_p : p.vertex2partition ) {
                std::pair<int64_t,int64_t> local_name = glo2loc[v_p.second];
                node_ids[local_name.first].push_back(v_p.second);
                partition[local_name.first].push_back(local_name.second);
            }

            for ( uint coarse_id = 0; coarse_id < p_coarse.num_partitions; ++coarse_id ) {
                std::vector<double> ideal_weights(local_partition_counter[coarse_id]);

                for ( uint i = 0; i < local_partition_counter[coarse_id]; ++i ) {
                    std::pair<int64_t,int64_t> local_name{ coarse_id, i };
                    ideal_weights[i] = p.ideal_weights[loc2glo[local_name]];
                }

                sub_partitions.emplace_back( local_partition_counter[coarse_id],
                                             g,
                                             node_ids[coarse_id],
                                             partition[coarse_id],
                                             ideal_weights);
            }
        }

        // add in moved vertices
        for ( int ID : moved_vertices ) {
            sub_partitions[ p_coarse.vertex2partition.at(ID) ].add_vertex(ID);
        }

        // refine local partitions
        for ( auto& part : sub_partitions ) {
            if ( part.num_partitions > 3 ) {
                std::cout << "Warning subpartition has more than 3 partitions. This will cause bad performance.\n"
                          << "     Consider lowering the coarsening factor to keep number of partitions in local\n"
                          << "     partitions low.\n";
            }

            refine_strict(g, part, 2);
            assert( part.is_balanced() );
        }

        //update big partition
        for ( uint coarse_id = 0; coarse_id < p_coarse.num_partitions; ++coarse_id ) {
            for ( auto& v_p : sub_partitions[coarse_id].vertex2partition ) {
                std::pair<int64_t,int64_t> local_name {coarse_id, v_p.second};
                p.vertex2partition[v_p.first] = loc2glo.at(local_name);
            }
        }
    }

    assert( p.is_balanced() );
}