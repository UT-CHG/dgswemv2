#include "../../partitioner/partition_type.hpp"
#include "../../partitioner/csrmat.hpp"

#include <numeric>

PartitionType coarsen(const CSRMat<>& g, const PartitionType& p, uint coarsening_factor);

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

    {
        std::vector<int> vertex_vec(8);
        std::vector<int64_t> p_vec(8);
        std::iota(vertex_vec.begin(), vertex_vec.end(), 0);
        std::iota(p_vec.begin(), p_vec.end(), 0);
        PartitionType p(8, g, vertex_vec, p_vec, std::vector<double>(8,1));

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

            double ideal_weight_sum{0};
            for ( auto& w : p2.ideal_weights ) {
                ideal_weight_sum += w;
            }
            if ( ideal_weight_sum != 8 ) {
                std::cerr << "Error: Ideal weights are not conserved across partitioning\n"
                          << "  Got " << ideal_weight_sum << " should be 8\n";
                error_found = true;
            }
        }
    }


    { // coarsen a subgraph
        std::vector<int> vertex_vec{0,1,2,4,5,6};
        std::vector<int64_t> p_vec(6);
        std::iota(p_vec.begin(), p_vec.end(), 0);
        PartitionType p(6, g, vertex_vec, p_vec, std::vector<double>(6,1));

        PartitionType p2 = coarsen(g,p,2);

        std::cout << "Partition for subgraph\n";
        for ( auto& v_p : p2.vertex2partition ) {
            std::cout << "  " << v_p.first << " : " << v_p.second << '\n';;
        }
        std::cout << '\n';


        if ( p2.num_partitions != 3 ) {
            std::cerr << "Error: Incorrect number of coarsening partitions have been generated\n"
                      << "  Got " << p2.num_partitions << " should be 3\n";
            error_found = true;
        }
        double ideal_weight_sum{0};
        for ( auto& w : p2.ideal_weights ) {
            ideal_weight_sum += w;
        }
        if ( ideal_weight_sum != 6 ) {
            std::cerr << "Error: Ideal weights are not conserved across partitioning\n"
                      << "  Got " << ideal_weight_sum << " should be 8\n";
            error_found = true;
        }
    }


    return error_found;
}