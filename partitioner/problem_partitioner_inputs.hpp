#ifndef PROBLEM_PARTITIONER_INPUTS_HPP
#define PROBLEM_PARTITIONER_INPUTS_HPP

#include <unordered_map>
#include <vector>

class ProblemPartitionerInputs {
  public:
    using WeightsType = std::unordered_map<int, std::vector<double>>;

    virtual WeightsType GetWeights() = 0;
    virtual void PartitionAuxiliaryFiles(){};

    virtual ~ProblemPartitionerInputs() = default;
};
#endif