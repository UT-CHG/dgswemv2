#ifndef IHDG_SWE_EDGE_DATA_STATE_HPP
#define IHDG_SWE_EDGE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeState {
    EdgeState() = default;
    EdgeState(const uint ndof, const uint ngp)
        : ze_hat(ndof),
          qx_hat(ndof),
          qy_hat(ndof),
          ze_hat_prev(ndof),
          qx_hat_prev(ndof),
          qy_hat_prev(ndof),
          ze_avg_at_gp(ngp),
          qx_avg_at_gp(ngp),
          qy_avg_at_gp(ngp),
          ze_hat_at_gp(ngp),
          qx_hat_at_gp(ngp),
          qy_hat_at_gp(ngp),
          h_hat_at_gp(ngp) {}

    std::vector<double> ze_hat;
    std::vector<double> qx_hat;
    std::vector<double> qy_hat;

    std::vector<double> ze_hat_prev;
    std::vector<double> qx_hat_prev;
    std::vector<double> qy_hat_prev;

    std::vector<double> ze_avg_at_gp;
    std::vector<double> qx_avg_at_gp;
    std::vector<double> qy_avg_at_gp;

    std::vector<double> ze_hat_at_gp;
    std::vector<double> qx_hat_at_gp;
    std::vector<double> qy_hat_at_gp;
    std::vector<double> h_hat_at_gp;
};
}
}

#endif
