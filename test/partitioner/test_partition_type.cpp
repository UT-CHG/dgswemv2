#include "../../partitioner/partition_type.hpp"

int main() {
    bool error_found{false};
    /*
     * Make a test graph that looks like
     *
     * 0--1--2--3
     * |  |  |  |
     * 5--6--7--8
     */

    std::unordered_map<int,double> nw;
    std::unordered_map<std::pair<int,int>,double> ew;

    for ( uint i = 0; i < 8; ++i ) {
        nw.insert({i,1.});
        if ( i < 7 ) {
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
    CSRMat<> meta_g = p.make_partition_graph(g);
    if ( !meta_g.get_node_wghts_map().count(0) || !meta_g.get_node_wghts_map().count(1) ) {
        std::cerr << "Can't find all of the super meshes\n";
        error_found = true;
    }

    if ( !meta_g.get_edge_wgts_map().count(std::make_pair(0,1)) ) {
        std::cerr << "Can't find edge between meshes\n";
        error_found = true;
    }

    return error_found;
}