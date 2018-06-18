#include "default_partitioner_inputs.hpp"

DefaultPartitionerInputs::DefaultPartitionerInputs(const MeshMetaData& mesh) {
    for (const auto& elt : mesh.elements) {
        weights.insert(std::make_pair(elt.first, std::vector<double>{1.}));
    }
}