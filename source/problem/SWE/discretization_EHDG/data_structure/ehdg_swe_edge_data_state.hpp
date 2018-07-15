#ifndef EHDG_SWE_EDGE_DATA_STATE_HPP
#define EHDG_SWE_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof) : q_hat(ndof), q_hat_prev(ndof) {}

    std::vector<StatVector<double, SWE::n_variables>> q_hat;
    std::vector<StatVector<double, SWE::n_variables>> q_hat_prev;
};
}
}

#endif
