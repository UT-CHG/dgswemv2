#ifndef GN_DATA_BOUNDARY_HPP
#define GN_DATA_BOUNDARY_HPP

namespace GN {
struct Boundary : SWE::Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : SWE::Boundary(ngp),
          dbath_hat_at_gp(GN::n_dimensions, ngp),
          w1_w1_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w1_w1_hat_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w2_w1_hat_kernel_at_gp(GN::n_dimensions, ngp),
          w1_hat_w1_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w1_hat_w2_kernel_at_gp(GN::n_dimensions, ngp) {}

    HybMatrix<double, GN::n_dimensions> dbath_hat_at_gp;

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

#endif
