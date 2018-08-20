#ifndef IHDG_GN_EDGE_DATA_STATE_HPP
#define IHDG_GN_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace GN {
namespace IHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : q_hat(GN::n_variables, ndof) {}

    /* swe containers */

    HybMatrix<double, GN::n_variables> q_hat;
};
}
}

#endif
