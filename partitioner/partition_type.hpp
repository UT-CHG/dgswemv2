#ifndef PARTITION_TYPE_HPP
#define PARTITION_TYPE_HPP

#include "csrmat.hpp"

#include <unordered_map>

struct PartitionType {
    uint num_partitions;
    std::unordered_map<int,int64_t> vertex2partition;

    PartitionType() = default;
    PartitionType(int64_t num_partitions, const CSRMat<>& g, const std::vector<int64_t>& partition);

    CSRMat<> make_partition_graph(const CSRMat<>& g);
};
#endif