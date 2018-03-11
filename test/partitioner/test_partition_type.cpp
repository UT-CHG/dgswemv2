#include "../../partitioner/partition_type.hpp"

int main() {
    bool error_found{false};
    /*
     * Make a test graph that looks like
     *
     * 0--1--2--3
     * |  |  |  |
     * 4--5--6--7
     */

    std::unordered_map<int,double> nw;
    std::unordered_map<std::pair<int,int>,double> ew;

    for ( uint i = 0; i < 8; ++i ) {
        nw.insert({i,1.});
        if ( i < 7 && i != 3) {
            std::pair<int,int> edge_name{i,i+1};
            ew.insert({edge_name,1});
        }

        if ( i < 4 ) {
            std::pair<int,int> edge_name{i,i+4};
            ew.insert({edge_name,1});
        }
    }

    CSRMat<> g(nw,ew);
    std::vector<std::function<double(int)>> cons;
    cons.push_back([](int id){ return 1;});

    std::vector<int64_t> p_vec = metis_part(g, 2, cons, 1.02);

    std::cout << "Partitioning for test graph:\n";
    for ( uint i = 0; i < 8;  ++i ) {
        std::cout << "  " << i << " : " << p_vec[i] << '\n';
    }
    {
        PartitionType p(2,g, p_vec);

        //check that the partitioning vector and partitioning map match up
        for ( uint indx = 0; indx < g.size(); ++indx ) {
            if ( p_vec[indx] != p.vertex2partition.at( g.node_ids()[indx]) ) {
                std::cerr << "Error in construction PartitionType\n"
                          << "Vertex " << g.node_ids()[indx] << " has been assigned to "
                          << p.vertex2partition.at( g.node_ids()[indx]) << " should be " << p_vec[indx] << '\n';
                error_found = true;
            }
        }

        //check that make_partition_graph is working properly
        CSRMat<> meta_g = p.make_partition_graph();
        if ( !meta_g.get_node_wghts_map().count(0) || !meta_g.get_node_wghts_map().count(1) ) {
            std::cerr << "Can't find all of the super meshes\n";
            error_found = true;
        }

        if ( !meta_g.get_edge_wgts_map().count(std::make_pair(0,1)) ) {
            std::cerr << "Can't find edge between meshes\n";
            error_found = true;
        }
    }


    { // test adding a vertex
        std::vector<int> vrtxs{0,2,4,5,6};
        std::vector<int64_t> part{0,1,0,0,1};

        PartitionType p(2,g,vrtxs,part,std::vector<double>(2,2));

        //add vertex one, which since this is greedy should get added to partition 0
        p.add_vertex(1);

        if ( p.vertex2partition[1] != 0 ) {
            std::cerr << "Error vertex 1 has been added to the incorrect partition\n";
            error_found = true;
        }
    }

    { // test is_balanced
        std::vector<std::vector<int>> vrtxs;
        std::vector<std::vector<int64_t>> parts;
        std::vector<std::vector<double>> ideal_weights;
        std::vector<bool> true_val;

        vrtxs.push_back({0,2,4,5,6});
        parts.push_back({0,1,0,0,1});
        ideal_weights.push_back({1,1});
        true_val.push_back(false);

        vrtxs.push_back({0,2,4,5,6});
        parts.push_back({0,1,0,0,1});
        ideal_weights.push_back({3,2});
        true_val.push_back(true);

        vrtxs.push_back({0,2,4,5,6});
        parts.push_back({0,1,0,0,1});
        ideal_weights.push_back({2,3});
        true_val.push_back(false);

        vrtxs.push_back({0,1,2,3,4,5,6,7});
        parts.push_back({0,0,1,1,0,0,1,1});
        ideal_weights.push_back({1,1});
        true_val.push_back(true);

        for ( uint case_=0; case_ < vrtxs.size(); ++case_ ) {
            PartitionType p(2,g,vrtxs[case_],parts[case_],ideal_weights[case_]);
            if ( p.is_balanced() != true_val[case_] ) {
                std::cerr << "Error in Partition::is_balanced() for case " << case_ << '\n';
                error_found = true;
            }
        }
    }

    return error_found;
}