#include "../../partitioner/partition_type.hpp"
#include "../../partitioner/csrmat.hpp"

void manually_refine2(const CSRMat<>& g, PartitionType& p, bool override=false);

bool test_manually_refine2(const CSRMat<>& g, PartitionType& p, const std::string& test_name) {
    bool error_found{false};

    manually_refine2(g,p);
    std::array<int,2> weights{0,0};
    for ( auto& v_p : p.vertex2partition ) {
        weights[v_p.second] += 1;
    }
    if ( weights[0] != weights[1] ) {
        std::cerr << "Error in test: " << test_name << '\n'
                  << "Weights in first partitioning are not equal.\n";
        error_found = true;
    }
    return error_found;
}

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

    { //steal one vertex
        {
            PartitionType p;

            p.num_partitions = 2;
            std::vector<int> p0{0,1,4}; // on partition 0
            std::vector<int> p1{2,3,5,6,7}; // on parititon 1
            for ( auto& v : p0 ) {
                p.vertex2partition.insert({v,0});
            }
            for ( auto& v : p1 ) {
                p.vertex2partition.insert({v,1});
            }
            error_found |= test_manually_refine2(g,p,"Steal 1 vertex");

            if ( p.vertex2partition[5] != 0 ) {
                std::cerr << "Error in test: Steal 1 vertex\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }

        }

        //flip partitions and try again
        {
            PartitionType p;

            p.num_partitions = 2;
            std::vector<int> p1{0,1,4}; // on partition 0
            std::vector<int> p0{2,3,5,6,7}; // on parititon 1
            for ( auto& v : p0 ) {
                p.vertex2partition.insert({v,0});
            }
            for ( auto& v : p1 ) {
                p.vertex2partition.insert({v,1});
            }
            error_found |= test_manually_refine2(g,p,"Steal 1 vertex flipped");

            if ( p.vertex2partition[5] != 1 ) {
                std::cerr << "Error in test: Steal 1 vertex flipped\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }

        }
    }

    { //steal two vertices
        {
            PartitionType p;

            p.num_partitions = 2;
            std::vector<int> p0{0,4}; // on partition 0
            std::vector<int> p1{1,2,3,5,6,7}; // on parititon 1
            for ( auto& v : p0 ) {
                p.vertex2partition.insert({v,0});
            }
            for ( auto& v : p1 ) {
                p.vertex2partition.insert({v,1});
            }
            error_found |= test_manually_refine2(g,p,"Steal 2 vertices");
            if ( p.vertex2partition[1] != 0 &&
                 p.vertex2partition[5] != 0 ) {
                std::cerr << "Error in test: Steal 2 vertices flipped\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }

        }

        //flip partitions and try again
        {
            PartitionType p;

            p.num_partitions = 2;
            std::vector<int> p1{0,4}; // on partition 0
            std::vector<int> p0{1,2,3,5,6,7}; // on parititon 1
            for ( auto& v : p0 ) {
                p.vertex2partition.insert({v,0});
            }
            for ( auto& v : p1 ) {
                p.vertex2partition.insert({v,1});
            }
            error_found |= test_manually_refine2(g,p,"Steal 2 vertices flipped");
            if ( p.vertex2partition[1] != 1 &&
                 p.vertex2partition[5] != 1 ) {
                std::cerr << "Error in test: Steal 2 vertices flipped\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }

        }
    }


    return error_found;
}