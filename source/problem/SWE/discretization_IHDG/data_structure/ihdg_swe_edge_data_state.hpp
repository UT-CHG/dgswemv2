#ifndef IHDG_SWE_EDGE_DATA_STATE_HPP
#define IHDG_SWE_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof)
        : ze_hat(ndof), qx_hat(ndof), qy_hat(ndof), ze_hat_prev(ndof), qx_hat_prev(ndof), qy_hat_prev(ndof) {}

    std::vector<double> ze_hat;
    std::vector<double> qx_hat;
    std::vector<double> qy_hat;

    std::vector<double> ze_hat_prev;
    std::vector<double> qx_hat_prev;
    std::vector<double> qy_hat_prev;
};
}
}

#endif
