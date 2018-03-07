#include "../partitioner/problem/default_partitioner_inputs.hpp"

#include "preprocessor/ADCIRC_reader/adcirc_format.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include <unordered_set>

std::vector<std::vector<MeshMetaData>> partition(const MeshMetaData& mesh_meta,
                                                 const std::unordered_map<int, std::vector<double>>& problem_weights,
                                                 const int num_partitions,
                                                 const int num_nodes,
                                                 const int ranks_per_locality,
                                                 const bool rank_balanced);

bool check_partition(const MeshMetaData& mesh, std::vector<std::vector<MeshMetaData>>& submeshes) {
    bool error_found{false};
    // 2 checks are performed in this unit test
    //  1.) Check that all the elements can be found in one of the submeshes
    {
        std::unordered_set<uint> elts;
        for (auto& e : mesh.elements) {
            elts.insert(e.first);
        }

        for (auto& loc : submeshes) {
            for (auto& sm : loc) {
                for (auto& e : sm.elements) {
                    elts.erase(e.first);
                }
            }
        }

        if (!elts.empty()) {
            error_found = true;
            std::cerr << "Error: Not all elements have been assigned to a submesh\n"
                      << "       Specifically elements: {";
            for (auto& e : elts) {
                std::cerr << " " << e;
            }
            std::cerr << " } have not been assigned\n";
        }
    }

    //  2.) Check that all the edges between submeshes are set correctly
    for (auto& loc : submeshes) {
        for (auto& sm : loc) {
            for (auto& e : sm.elements) {
                for (uint k = 0; k < e.second.neighbor_ID.size(); ++k) {
                    if (e.second.boundary_type[k] == INTERNAL && !sm.elements.count(e.second.neighbor_ID[k])) {
                        error_found = true;
                        std::cerr << "Element " << e.first << " cannot find internal neighbor of id "
                                  << e.second.neighbor_ID[k] << '\n';
                    } else if (e.second.boundary_type[k] == DISTRIBUTED && sm.elements.count(e.second.neighbor_ID[k])) {
                        error_found = true;
                        std::cerr << "Element " << e.first << " found distributed neighbor on its own submeh\n";
                    }
                }
            }
        }
    }

    return error_found;
}

// This unit test partitions a mesh and then stitches it back together again
int main(int argc, char** argv) {

    bool error_found{false};

    AdcircFormat mesh1(argv[1]);
    MeshMetaData meshA(mesh1);

    std::unordered_map<int,std::vector<double>> weights;
    for ( const auto& elt : meshA.elements ) {
        weights.insert(std::make_pair(elt.first, std::vector<double>{1.}));
    }

    std::cout << "Checking for 1 locality...\n";
    {
        std::vector<std::vector<MeshMetaData>> submeshes = partition(meshA, weights, 4, 1, 1, false);
        error_found = error_found || check_partition(meshA, submeshes);
    }
    std::cout << "...done checking for 1 locality\n\n";
    std::cout << "Checking for multiple localities with one rank per locality...\n";
    {
        std::vector<std::vector<MeshMetaData>> submeshes = partition(meshA, weights, 4, 2, 1, false);
        error_found = error_found || check_partition(meshA, submeshes);
    }
    std::cout << "...done checking for multiple localities\n";
    std::cout << "Checking for multiple localities with multiple ranks per locality...\n";
    {
        std::vector<std::vector<MeshMetaData>> submeshes = partition(meshA, weights, 8, 2, 2, false);
        error_found = error_found || check_partition(meshA, submeshes);
    }
    std::cout << "...done checking for multiple localities\n";
    std::cout << "Checking for with rank balancing turned on...\n";
    {
        std::vector<std::vector<MeshMetaData>> submeshes = partition(meshA, weights, 4, 1, 1, true);
        error_found = error_found || check_partition(meshA, submeshes);
    }
    std::cout << "...done checking for rank balancing turned on\n\n";

    if (error_found) {
        return 1;
    }

    return 0;
}
