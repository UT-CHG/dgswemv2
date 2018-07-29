#ifndef EHDG_GN_EDGE_DATA_STATE_HPP
#define EHDG_GN_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : q_hat(GN::n_variables, ndof) {}

    HybMatrix<double, GN::n_variables> q_hat;
};
}
}

#endif
