#ifndef EHDG_GN_DATA_BOUNDARY_HPP
#define EHDG_GN_DATA_BOUNDARY_HPP

namespace GN {
namespace EHDG {
struct Boundary : SWE::EHDG::Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : SWE::EHDG::Boundary(ngp),
          ze_hat_at_gp(ngp),
          bath_hat_at_gp(ngp),
          dbath_hat_at_gp(GN::n_dimensions, ngp),
          ddbath_hat_at_gp(GN::n_ddbath_terms, ngp),
          u_hat_at_gp(GN::n_dimensions, ngp),
          du_hat_at_gp(GN::n_du_terms, ngp),
          w1_w1_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w1_w1_hat_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w2_w1_hat_kernel_at_gp(GN::n_dimensions, ngp),
          w1_hat_w1_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w1_hat_w2_kernel_at_gp(GN::n_dimensions, ngp) {}

    DynRowVector<double> ze_hat_at_gp;

    DynRowVector<double> bath_hat_at_gp;
    HybMatrix<double, GN::n_dimensions> dbath_hat_at_gp;
    HybMatrix<double, GN::n_ddbath_terms> ddbath_hat_at_gp;

    HybMatrix<double, GN::n_dimensions> u_hat_at_gp;
    HybMatrix<double, GN::n_du_terms> du_hat_at_gp;

    HybMatrix<double, GN::n_dimensions * GN::n_dimensions> w1_w1_kernel_at_gp;

    HybMatrix<double, GN::n_dimensions * GN::n_dimensions> w1_w1_hat_kernel_at_gp;
    HybMatrix<double, GN::n_dimensions> w2_w1_hat_kernel_at_gp;

    HybMatrix<double, GN::n_dimensions * GN::n_dimensions> w1_hat_w1_kernel_at_gp;
    HybMatrix<double, GN::n_dimensions> w1_hat_w2_kernel_at_gp;

    DynMatrix<double> w1_w1_hat;
    DynMatrix<double> w2_w1_hat;

    DynMatrix<double> w1_hat_w1;
    DynMatrix<double> w1_hat_w2;

    std::vector<uint> dc_global_dof_indx;
};
}
}

#endif
