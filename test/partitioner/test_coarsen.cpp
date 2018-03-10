#include "../../partitioner/partition_type.hpp"
#include "../../partitioner/csrmat.hpp"

#include <numeric>

PartitionType coarsen(const CSRMat<>& g, PartitionType& p, uint coarsening_factor);

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
    std::vector<int64_t> p_vec(8);
    std::iota(p_vec.begin(), p_vec.end(), 0);
    PartitionType p(8, g, p_vec);

    std::vector<uint> cf{2,3,4};
    std::vector<uint> true_num_partitions{4,3,2};

    for ( uint i = 0; i < cf.size(); ++i ) {
        PartitionType p2 = coarsen(g, p, cf[i]);

        std::cout << "Partition for cf = " << cf[i] << '\n';
        for ( uint j = 0; j < 8; ++j ) {
            std::cout << "  " << j << " : " << p2.vertex2partition.at(j) << '\n';;
        }
        std::cout << '\n';

        if ( p2.num_partitions != true_num_partitions[i] ) {
            std::cerr << "Error: Incorrect number of coarsening partitions have been generated\n"
                      << "  Got " << p2.num_partitions << " should be " << true_num_partitions[i] << '\n';
            error_found = true;
        }
    }

    return error_found;
}