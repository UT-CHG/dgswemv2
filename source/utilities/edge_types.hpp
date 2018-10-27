#ifndef EDGE_TYPES_HPP
#define EDGE_TYPES_HPP

constexpr uchar INTERNAL{1 << 7};            // 128
constexpr uchar DISTRIBUTED_OFFSET{1 << 6};  // 64
constexpr uchar DISTRIBUTED{INTERNAL + DISTRIBUTED_OFFSET};

constexpr bool is_internal(uchar edge) {
    return edge >= INTERNAL && edge < DISTRIBUTED;
}

constexpr bool is_distributed(uchar edge) {
    return edge >= DISTRIBUTED;
}

constexpr bool is_boundary(uchar edge) {
    return edge < INTERNAL;
}

constexpr uchar distributed(uchar internal_edge) {
    return internal_edge + DISTRIBUTED_OFFSET;
}

#endif