#ifndef SWE_PARTITIONER_INPUTS_HPP
#define SWE_PARTITIONER_INPUTS_HPP

#include "../problem_partitioner_inputs.hpp"

#include "preprocessor/mesh_metadata.hpp"
#include "problem/SWE/problem_input/swe_inputs.hpp"

namespace SWE {
class PartitionerInputs final : public ProblemPartitionerInputs {
  public:
    PartitionerInputs(const MeshMetaData& mesh, Inputs inputs);

    inline WeightsType GetWeights() { return weights; }

    // To be added at a later date
    // void PartitionAuxiliaryFiles();

  private:
    WeightsType weights;
};
}
#endif