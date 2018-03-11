#ifndef PARTITION_TYPE_HPP
#define PARTITION_TYPE_HPP

#include "csrmat.hpp"

#include <unordered_map>

struct PartitionType {
    uint num_partitions;
    const CSRMat<>& g_ref;
    std::unordered_map<int,int64_t> vertex2partition;
    std::vector<double> ideal_weights;

    PartitionType() = default;

    PartitionType(int64_t num_partitions,
                  const CSRMat<>& g,
                  const std::vector<int64_t>& partition);

    PartitionType(int64_t num_partitions,
                  const CSRMat<>& g,
                  const std::vector<int64_t>& partition,
                  const std::vector<double>& ideal_weights);

    PartitionType(int64_t num_partitions,
                  const CSRMat<>& g,
                  const std::vector<int>& node_ids,
                  const std::vector<int64_t>& partition,
                  const std::vector<double>& ideal_weights);

    //Constructor to generate coarse partition of p based on coarse partition
    PartitionType(int64_t num_partitions,
                  const CSRMat<>& g,
                  const PartitionType& p,
                  const std::vector<int64_t>& coarse_partition);

    CSRMat<> make_partition_graph() const;

    //needs test//implementation
    bool is_balanced();

    //needs test
    void add_vertex(int ID);
};
#endif