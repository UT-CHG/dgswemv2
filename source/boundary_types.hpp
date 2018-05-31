#ifndef BOUNDARY_TYPES_HPP
#define BOUNDARY_TYPES_HPP

constexpr uchar INTERNAL{1<<7}; //128
constexpr uchar DISTRIBUTED_OFFSET{1<<6}; //64
constexpr uchar DISTRIBUTED{INTERNAL+DISTRIBUTED_OFFSET};

constexpr bool IsInternal(uchar boundary) {
    return boundary >= INTERNAL && boundary < DISTRIBUTED;
}

constexpr bool IsDistributed(uchar boundary) {
    return boundary >= DISTRIBUTED;
}

constexpr bool IsBoundary(uchar boundary) {
    return boundary < INTERNAL;
}


template <uchar Boundary>
struct Distributed {
    static constexpr uchar Type() {
        static_assert(IsInternal(Boundary),
            "Error: Distribute<Boundary>::Type() must be an internal boundary, which ranges from INTERNAL to DISTRIBUTED");
        return Boundary + DISTRIBUTED_OFFSET;
    }
};
#endif