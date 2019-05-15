#ifndef DEFAULT_PARTITIONER_INPUTS_HPP
#define DEFAULT_PARTITIONER_INPUTS_HPP

#include "../problem_partitioner_inputs.hpp"

#include "preprocessor/mesh_metadata.hpp"

class DefaultPartitionerInputs final : public ProblemPartitionerInputs {
  public:
    DefaultPartitionerInputs(const MeshMetaData& mesh);

    WeightsType GetWeights() final { return weights; }

  private:
    WeightsType weights;
};
#endif