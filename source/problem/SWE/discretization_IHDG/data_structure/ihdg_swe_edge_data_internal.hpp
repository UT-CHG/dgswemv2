#ifndef IHDG_SWE_EDGE_DATA_INTERNAL_HPP
#define IHDG_SWE_EDGE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : ze_avg_at_gp(ngp),
          qx_avg_at_gp(ngp),
          qy_avg_at_gp(ngp),
          ze_hat_at_gp(ngp),
          qx_hat_at_gp(ngp),
          qy_hat_at_gp(ngp),
          h_hat_at_gp(ngp) {}

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
