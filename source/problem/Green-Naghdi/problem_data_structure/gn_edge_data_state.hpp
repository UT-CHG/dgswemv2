#ifndef GN_EDGE_DATA_STATE_HPP
#define GN_EDGE_DATA_STATE_HPP

namespace GN {
struct EdgeState : SWE::EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : SWE::EdgeState(ndof) {}
};
}

#endif
