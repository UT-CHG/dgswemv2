#ifndef IHDG_SWE_EDGE_DATA_STATE_HPP
#define IHDG_SWE_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : q_hat(SWE::n_variables, ndof), q_hat_prev(SWE::n_variables, ndof) {}

    DynMatrix<double> q_hat;
    DynMatrix<double> q_hat_prev;
};
}
}

#endif
