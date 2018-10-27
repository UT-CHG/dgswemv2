#ifndef IHDG_SWE_EDGE_DATA_STATE_HPP
#define IHDG_SWE_EDGE_DATA_STATE_HPP

namespace SWE {
namespace IHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : q_hat(SWE::n_variables, ndof) {}

    HybMatrix<double, SWE::n_variables> q_hat;
};
}
}

#endif
