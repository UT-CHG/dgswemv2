#ifndef PARTITION_TYPE_HPP
#define PARTITION_TYPE_HPP

#include <unordered_map>

struct PartitionType {
    uint num_partitions;
    std::unordered_map<int,int64_t> vertex2partition;
};
#endif