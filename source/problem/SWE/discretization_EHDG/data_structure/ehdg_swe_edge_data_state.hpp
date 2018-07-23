#ifndef EHDG_SWE_EDGE_DATA_STATE_HPP
#define EHDG_SWE_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : q_hat(SWE::n_variables, ndof), q_hat_prev(SWE::n_variables, ndof) {}

    HybMatrix<double, SWE::n_variables> q_hat;
    HybMatrix<double, SWE::n_variables> q_hat_prev;
};
}
}

#endif
