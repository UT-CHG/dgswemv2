#include "partition_type.hpp"
#include "csrmat.hpp"
#include "steal_queue.hpp"

#include "utilities/almost_equal.hpp"

#include <deque>
#include <memory>

uint spaces = 0;

std::string spacing() {
    std::string white_space;
    for ( uint i = 0; i < spaces; ++i ) {
        white_space += " ";
    }
    return white_space;
}

void manually_refine2(const CSRMat<>& g, PartitionType&p,
                      const std::vector<double>& custom_weights = std::vector<double>()) {
    spaces += 2;
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

    std::cout << spacing() << "high: " << high << " low: " << low << '\n';

    std::cout << spacing() << "High elements:";
    for( auto& elt : elements[high] ) {
        std::cout << ' ' << elt;
    }
    std::cout << '\n';

    StealQ steal_q(g, p, elements[high]);

    std::cout << spacing() << "num_tiles_to_be_moved: " << num_tiles_to_be_moved << '\n';

    while ( num_tiles_to_be_moved > 0 ) {
        int id = steal_q.steal_one();
        std::cout <<spacing() << "stealing id: " << id << '\n';
        assert( p.vertex2partition[id] == high );
        p.vertex2partition[id] = low;
        --num_tiles_to_be_moved;
    }
    spaces -= 2;
}
/*
void manually_refine3() {

}

*/
PartitionType coarsen(const CSRMat<>& g, const PartitionType& p, uint coarsening_factor) {
  //easy case
  if ( coarsening_factor == 1 ) {
    return p;
  }

  //if exact these are the number of partitions
    int64_t num_partitions = p.vertex2partition.size()/coarsening_factor;

    //Otherwise, add an extra partition
    if ( num_partitions*coarsening_factor != p.num_partitions ) {
        num_partitions = static_cast<int64_t>( std::ceil( static_cast<double>(p.num_partitions)/coarsening_factor));
    }
    assert( num_partitions <= g.size() );

    std::cout << spacing() << "num_partitions: " << num_partitions << '\n';

    CSRMat<> g_super = p.make_partition_graph();
    std::cout << spacing() << "g_super.size() = " << g_super.size() << '\n';
    
    std::vector<std::function<double(int)>> cons;
    cons.push_back([&g_super](int id) { return g_super.node_weight(id); });

    std::vector<int64_t> coarse_partition = metis_part(g_super, num_partitions, cons, 1.02);

    return PartitionType(num_partitions,
                         g,
                         p,
                         g_super.node_ids(),
                         coarse_partition);
}

void refine_strict(const CSRMat<>& g, PartitionType& p, uint coarsening_factor) {
  assert( std::addressof(g) == std::addressof(p.g_ref) );

    if ( coarsening_factor <= 1 ) {
        throw std::logic_error("Error in refine_strict;"
                               "Coarsening factor needs to be greater than 1\n");
    }

    //Nothing to be done
    if ( p.is_balanced() ) {
        std::cout << spacing() << "Partition is balanced; Nothing to be done\n\n";
        return;
    }

    //Take care of base cases
    if ( p.num_partitions == 1 ) {
        std::cout << spacing() << "One partition; Trivially balanced\n\n";
        return;
    } else if ( p.num_partitions == 2 ) {
        std::cout << spacing() << "Two partitions; manually refining\n\n";
        manually_refine2(g,p);
        return;
/*    } else if ( p.num_partitions ==3 ) {
        manually_refine3(g,p);
        return;*/
    }


    std::cout << spacing() << "Coarsening partition..." << coarsening_factor <<
 "\n";
    PartitionType p_coarse = coarsen(g, p, coarsening_factor);
    //std::cout << "p_coarse summary:\n" << p_coarse;
    std::unordered_set<int> moved_vertices;


    { //save copy of old vertex placements to track movement
        PartitionType p_coarse_old = p_coarse;
        std::cout << spacing() << "Entering refine strict...\n";
        spaces += 2;
        refine_strict(g, p_coarse, 2);
        spaces -= 2;
        std::cout << spacing() << "Leaving refine strict...\n";

        for ( auto& v_p : p_coarse.vertex2partition ) {
            if ( v_p.second != p_coarse_old.vertex2partition.at(v_p.first) ) {
                moved_vertices.insert(v_p.first);
                //moved vertices are temporarily not assinged to any partition
                p.vertex2partition.erase(v_p.first);
            }
        }

        std::cout << spacing() << "moved vertices:";
        for ( auto& v : moved_vertices ) {
            std::cout << " " << v;
        }
        std::cout << '\n';
    }

    { //make sub partitioned subgraphs
        std::vector<std::pair<int64_t,int64_t>> glo2loc(p.num_partitions,std::make_pair(-1,-1));
        std::unordered_map<std::pair<int64_t,int64_t>, int64_t> loc2glo;
        std::vector<std::size_t> local_partition_counter(p_coarse.num_partitions,0);

        for ( auto& v_p : p.vertex2partition ) {
            if ( glo2loc[v_p.second].first == -1 ) { //hasn't been assigned yet
                int64_t coarse_partition = p_coarse.vertex2partition[v_p.first];
                int64_t local_partition = local_partition_counter[coarse_partition]++;
                std::pair<int64_t,int64_t> local_name{coarse_partition, local_partition};

                glo2loc[v_p.second] = local_name;
                loc2glo[local_name] = v_p.second;
            }
        }

        /*std::cout << spacing() << "glo2loc:\n";
        for ( uint i = 0; i < glo2loc.size(); ++i ) {
            std::cout << spacing() << "  " << i << " : " << glo2loc[i].first << ", " << glo2loc[i].second << '\n';
        }

        std::cout << spacing() << "loc2glo:\n";
        for ( auto& l_p : loc2glo ) {
            std::cout << spacing() << "  " << l_p.first.first << ", " << l_p.first.second << " : " << l_p.second << '\n';
            }*/

        //build the local sub_partitions
        std::vector<PartitionType> sub_partitions;
        {
            std::vector<std::vector<int>> node_ids(p_coarse.num_partitions);
            std::vector<std::vector<int64_t>> partition(p_coarse.num_partitions);
            for ( auto& v_p : p.vertex2partition ) {
                std::pair<int64_t,int64_t> local_name = glo2loc.at(v_p.second);
                node_ids[local_name.first].push_back(v_p.first);
                partition[local_name.first].push_back(local_name.second);
            }

            sub_partitions.reserve(p_coarse.num_partitions);
            for ( uint coarse_id = 0; coarse_id < p_coarse.num_partitions; ++coarse_id ) {
                std::vector<double> ideal_weights(local_partition_counter[coarse_id]);

                for ( uint i = 0; i < local_partition_counter[coarse_id]; ++i ) {
                    std::pair<int64_t,int64_t> local_name{ coarse_id, i };
                    ideal_weights[i] = p.ideal_weights.at(loc2glo.at(local_name));
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

        /*for ( uint sp = 0; sp < sub_partitions.size(); ++sp ) {
          std::cout << spacing()<< "Subpartition: " << sp << '\n' << sub_partitions[sp];
          }*/

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

void refine_strict_tracked(const CSRMat<>& g, std::vector<PartitionType>& ps, const std::vector<uint>& coarsening_factors) {
    int recursion_level = ps.size() - 1;
    std::cout << spacing() << "recursion_level: " << recursion_level << std::endl;

    std::cout << spacing() << "current partition is_balanced() = " << std::boolalpha << ps[recursion_level].is_balanced() << '\n';
    std::cout << spacing() << "Coarsening partition\n";
    ps.emplace_back( coarsen(g, ps[recursion_level], coarsening_factors[recursion_level]) );

    //std::cout << "p_coarse summary:\n" << p_coarse;
    std::unordered_set<int> moved_vertices;


    { //save copy of old vertex placements to track movement
        PartitionType p_coarse_old = ps[recursion_level + 1];
        if ( recursion_level < coarsening_factors.size() -1 ) {
          std::cout << spacing() << "Entering refine strict tracked...\n";
          spaces += 2;
          refine_strict_tracked(g, ps, coarsening_factors);
          spaces -= 2;
          std::cout << spacing() << "Leaving refine strict tracked...\n";

        } else {
          std::cout << spacing() << "Entering refine strict...\n";
          spaces += 2;
          refine_strict(g, ps[recursion_level + 1], 2);
          spaces -= 2;
          std::cout << spacing() << "Leaving refine strict...\n";
        }

        for ( auto& v_p : ps[recursion_level + 1].vertex2partition ) {
            if ( v_p.second != p_coarse_old.vertex2partition.at(v_p.first) ) {
                moved_vertices.insert(v_p.first);
                //moved vertices are temporarily not assinged to any partition
                ps[recursion_level].vertex2partition.erase(v_p.first);
            }
        }

        std::cout << spacing() << "moved vertices:";
        for ( auto& v : moved_vertices ) {
            std::cout << " " << v;
        }
        std::cout << '\n';
    }

    { //make sub partitioned subgraphs
        std::vector<std::pair<int64_t,int64_t>> glo2loc(ps[recursion_level].num_partitions,std::make_pair(-1,-1));
        std::unordered_map<std::pair<int64_t,int64_t>, int64_t> loc2glo;
        std::vector<std::size_t> local_partition_counter(ps[recursion_level + 1].num_partitions,0);

        for ( auto& v_p : ps[recursion_level].vertex2partition ) {
            if ( glo2loc[v_p.second].first == -1 ) { //hasn't been assigned yet
                int64_t coarse_partition = ps[recursion_level + 1].vertex2partition[v_p.first];
                int64_t local_partition = local_partition_counter[coarse_partition]++;
                std::pair<int64_t,int64_t> local_name{coarse_partition, local_partition};

                glo2loc[v_p.second] = local_name;
                loc2glo[local_name] = v_p.second;
            }
        }

        std::cout << spacing() << "glo2loc:\n";
        for ( uint i = 0; i < glo2loc.size(); ++i ) {
            std::cout << spacing() << "  " << i << " : " << glo2loc[i].first << ", " << glo2loc[i].second << '\n';
        }

        std::cout << spacing() << "loc2glo:\n";
        for ( auto& l_p : loc2glo ) {
            std::cout << spacing() << "  " << l_p.first.first << ", " << l_p.first.second << " : " << l_p.second << '\n';
        }

        //build the local sub_partitions
        std::vector<PartitionType> sub_partitions;
        {
            std::vector<std::vector<int>> node_ids(ps[recursion_level+1].num_partitions);
            std::vector<std::vector<int64_t>> partition(ps[recursion_level+1].num_partitions);
            for ( auto& v_p : ps[recursion_level].vertex2partition ) {
                std::pair<int64_t,int64_t> local_name = glo2loc.at(v_p.second);
                node_ids[local_name.first].push_back(v_p.first);
                partition[local_name.first].push_back(local_name.second);
            }

            sub_partitions.reserve(ps[recursion_level+1].num_partitions);
            for ( uint coarse_id = 0; coarse_id < ps[recursion_level+1].num_partitions; ++coarse_id ) {
                std::vector<double> ideal_weights(local_partition_counter[coarse_id]);

                for ( uint i = 0; i < local_partition_counter[coarse_id]; ++i ) {
                    std::pair<int64_t,int64_t> local_name{ coarse_id, i };
                    ideal_weights[i] = ps[recursion_level].ideal_weights.at(loc2glo.at(local_name));
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
            sub_partitions[ ps[recursion_level+1].vertex2partition.at(ID) ].add_vertex(ID);
        }

        /*for ( uint sp = 0; sp < sub_partitions.size(); ++sp ) {
            std::cout << "Subpartition: " << sp << '\n' << sub_partitions[sp];
            }*/

        // refine local partitions
        for ( auto& part : sub_partitions ) {
            refine_strict(g, part, 2);
            assert( part.is_balanced() );
        }

        //update big partition
        for ( uint coarse_id = 0; coarse_id < ps[recursion_level+1].num_partitions; ++coarse_id ) {
            for ( auto& v_p : sub_partitions[coarse_id].vertex2partition ) {
                std::pair<int64_t,int64_t> local_name {coarse_id, v_p.second};
                ps[recursion_level].vertex2partition[v_p.first] = loc2glo.at(local_name);
            }
        }
    }

    assert( ps[recursion_level].is_balanced() );
}
