#ifndef GN_EDGE_DATA_INTERNAL_HPP
#define GN_EDGE_DATA_INTERNAL_HPP

namespace GN {
struct EdgeInternal : SWE::EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : SWE::EdgeInternal(ngp), w1_hat_w1_hat_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp) {}

    HybMatrix<double, GN::n_dimensions * GN::n_dimensions> w1_hat_w1_hat_kernel_at_gp;

    DynMatrix<double> w1_hat_w1_hat;
    DynVector<double> w1_hat_rhs;

    DynVector<double> w1_hat_w1_hat_flat;
    std::vector<DynVector<double>> w1_hat_w1_hat_con_flat;

    std::vector<uint> dc_global_dof_indx;
    uint dc_sol_offset;
};
}

#endif
