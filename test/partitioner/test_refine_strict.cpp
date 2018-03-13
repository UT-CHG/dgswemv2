#include "../../partitioner/partition_type.hpp"
#include "../../partitioner/csrmat.hpp"

void refine_strict(const CSRMat<>& g, PartitionType& p, uint coarsening_factor);

CSRMat<> make_test_graph(uint n) {
    /*
     * Make a test graph that looks like
     *
     * 0---1-----2---...--(n-1)
     * |   |     |          |
     * n-(n+1)-(n+2)--...--(2*n-1)
     */

    std::unordered_map<int,double> nw;
    std::unordered_map<std::pair<int,int>,double> ew;

    for ( uint i = 0; i < 2*n; ++i ) {
        nw.insert({i,1.});
        if ( (i < (2*n-1) ) && (i != n-1)) {
            std::pair<int,int> edge_name{i,i+1};
            ew.insert({edge_name,1});
        }

        if ( i < n ) {
            std::pair<int,int> edge_name{i,i+n};
            ew.insert({edge_name,1});
        }
    }

    return CSRMat<>(nw,ew);
}

int main(int argc, char** argv) {
    CSRMat<>g = make_test_graph(6);

    std::vector<int64_t> p_vec{0,0,1,1,2,2,0,1,1,2,2,2};

    std::cout << "g.size() = " << g.size() << '\n';

    g.csr_info();

    std::cout << "p_vec.size() = " << p_vec.size() << '\n';

    PartitionType p(3, g, p_vec);
    std::cout << "initial partition\n" << p << '\n';


    refine_strict(g,p,2);

    std::cout << "Post strict refine partition\n";
    for ( uint vrtx = 0; vrtx < 12; ++vrtx ) {
        std::cout << "  " << vrtx << " : " << p.vertex2partition[vrtx] << '\n';
    }

    return !p.is_balanced();
}