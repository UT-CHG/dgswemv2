#ifndef EHDG_GN_EDGE_DATA_STATE_HPP
#define EHDG_GN_EDGE_DATA_STATE_HPP

namespace GN {
namespace EHDG {
struct EdgeState : SWE_SIM::EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : SWE_SIM::EdgeState(ndof) {}
};
}
}

#endif
