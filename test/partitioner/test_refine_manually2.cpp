#include "../../partitioner/partition_type.hpp"
#include "../../partitioner/csrmat.hpp"

#include "utilities/ignore.hpp"

void manually_refine2(const CSRMat<>& g, PartitionType& p,
                      const std::vector<double>& custom_weights = std::vector<double>());

bool test_manually_refine2(const CSRMat<>& g, PartitionType& p, const std::string& test_name) {
    bool error_found{false};

    manually_refine2(g,p);
    if ( !p.is_balanced() ) {
        std::cerr << "Error in test: " << test_name << '\n'
                  << "Weights in first partitioning are not equal.\n";
        error_found = true;
    }
    return error_found;
}

void print_partition_type(const PartitionType& p) {
    for ( const auto& v_p : p.vertex2partition ) {
        std::cout << v_p.first << " : " << v_p.second << '\n';
    }
    std::cout << '\n';
}

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
        if ( i < 7 && i !=3 ) {
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
            std::vector<int64_t> p_vec{0,0,1,1,0,1,1,1};
            PartitionType p(2,g,p_vec);;

            error_found |= test_manually_refine2(g,p,"Steal 1 vertex");

            if ( p.vertex2partition[5] != 0 ) {
                std::cerr << "Error in test: Steal 1 vertex\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }
        }

        //flip partitions and try again
        {
            std::vector<int64_t> p_vec{1,1,0,0,1,0,0,0};
            PartitionType p(2,g,p_vec);

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
            std::vector<int64_t> p_vec{0,1,1,1,0,1,1,1};
            PartitionType p(2,g,p_vec);

            error_found |= test_manually_refine2(g,p,"Steal 2 vertices");
            if ( p.vertex2partition[1] != 0 &&
                 p.vertex2partition[5] != 0 ) {
                std::cerr << "Error in test: Steal 2 vertices flipped\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }
        }

        { //flip partitions and try again
            std::vector<int64_t> p_vec{1,0,0,0,1,0,0,0};
            PartitionType p(2,g,p_vec);

            error_found |= test_manually_refine2(g,p,"Steal 2 vertices flipped");
            if ( p.vertex2partition[1] != 1 &&
                 p.vertex2partition[5] != 1 ) {
                std::cerr << "Error in test: Steal 2 vertices flipped\n"
                          << "Didn't send over the correct vertex.\n";
                error_found = true;
            }

        }
    }

    { // refine on a subgraph
        std::vector<int> vrtxs{1,2,5,6};
        std::vector<int64_t> p_vec{1,0,0,0};
        std::vector<double> iw{1,1};
        PartitionType p(2,g,vrtxs,p_vec,iw);

        error_found |= test_manually_refine2(g,p,"Steal on subgraph");
        if ( p.vertex2partition[2] != 1 &&
             p.vertex2partition[5] != 1 ) {
            print_partition_type(p);

            std::cerr << "Error in test: Steal on subgraph\n"
                      << "Didn't send over the correct vertex.\n";
            error_found = true;
        }
    }

    { // steal with non-equal ideal weights
        std::vector<int> vrtxs{0,1,2,3,4,5,6,7};
        std::vector<int64_t> p_vec{0,0,0,1,0,0,1,1};
        std::vector<double> iw{3,1};
        PartitionType p(2,g,vrtxs,p_vec,iw);

        error_found |= test_manually_refine2(g,p,"Steal with non-trivial ideal weights");
        if ( p.vertex2partition[6] != 0 ) {
            print_partition_type(p);

            std::cerr << "Error in test: Steal with non-trivial ideal weights\n"
                      << "Didn't send over the correct vertex.\n";
            error_found = true;
        }
    }

    { //try an impossible steal
        bool local_error{true};
        try {
            std::vector<int> vrtxs{1,2,3,5,6};
            std::vector<int64_t> p_vec{1,0,0,0,0};
            std::vector<double> iw{1,1};
            PartitionType p(2,g,vrtxs,p_vec,iw);

            bool dummy = test_manually_refine2(g,p,"Steal on subgraph");
            Utilities::ignore(dummy);
        } catch (...) {
            local_error = false;
        }
        if (local_error) {
            std::cerr << "Error: manually_refine2 attempted an impossible steal\n"
                      << "       and did not throw an error.\n";
            error_found = true;
        }
    }

    { //steal with custom weights
        std::vector<int> vrtxs{1,2,3,5,6,7};
        std::vector<int64_t> p_vec{1,1,0,1,0,0};
        std::vector<double> iw{1,1};
        PartitionType p(2,g,vrtxs,p_vec,iw);

        manually_refine2(g,p, std::vector<double>{1,2});

        if ( p.vertex2partition[6] != 1 ) {
            print_partition_type(p);

            std::cerr << "Error in test: Steal with custom weights\n"
                      << "Didn't send over the correct vertex.\n";
            error_found = true;
        }
    }

    return error_found;
}