#ifndef EHDG_GN_EDGE_DATA_INTERNAL_HPP
#define EHDG_GN_EDGE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct EdgeInternal : SWE::EHDG::EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : SWE::EHDG::EdgeInternal(ngp), w1_hat_w1_hat_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp) {}

    HybMatrix<double, GN::n_dimensions * GN::n_dimensions> w1_hat_w1_hat_kernel_at_gp;

    DynMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;

    std::vector<uint> dc_global_dof_indx;
    uint dc_sol_offset;
};
}
}

#endif
