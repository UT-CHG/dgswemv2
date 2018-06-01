#ifndef EDGE_TYPES_HPP
#define EDGE_TYPES_HPP

using uchar = unsigned char;

constexpr uchar INTERNAL{1<<7}; //128
constexpr uchar DISTRIBUTED_OFFSET{1<<6}; //64
constexpr uchar DISTRIBUTED{INTERNAL+DISTRIBUTED_OFFSET};

constexpr bool IsInternal(uchar edge) {
    return edge >= INTERNAL && edge < DISTRIBUTED;
}

constexpr bool IsDistributed(uchar edge) {
    return edge >= DISTRIBUTED;
}

constexpr bool IsBoundary(uchar edge) {
    return edge < INTERNAL;
}


template <uchar EdgeType>
struct Distributed {
    static constexpr uchar Type() {
        static_assert(IsInternal(EdgeType),
            "Error: Distribute<EdgeTYpe>::Type() must be an internal boundary, which ranges from INTERNAL to DISTRIBUTED");
        return EdgeType + DISTRIBUTED_OFFSET;
    }
};
#endif