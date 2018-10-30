#ifndef EHDG_GN_DATA_INTERNAL_HPP
#define EHDG_GN_DATA_INTERNAL_HPP

namespace GN {
namespace EHDG {
struct Internal : SWE_SIM::Internal {
    Internal() = default;
    Internal(const uint ngp)
        : SWE_SIM::Internal(ngp),
          u_at_gp(GN::n_dimensions, ngp),
          du_at_gp(GN::n_du_terms, ngp),
          ddu_at_gp(GN::n_ddu_terms, ngp),
          ddbath_at_gp(GN::n_ddbath_terms, ngp),
          dddbath_at_gp(GN::n_dddbath_terms, ngp),
          w1_w1_kernel_at_gp(GN::n_dimensions * GN::n_dimensions, ngp),
          w1_w2_kernel_at_gp(GN::n_dimensions, ngp),
          w1_rhs_kernel_at_gp(GN::n_dimensions, ngp),
          w2_w1_kernel_at_gp(GN::n_dimensions, ngp),
          w2_w2_kernel_at_gp(ngp) {}

    HybMatrix<double, GN::n_dimensions> u_at_gp;
    HybMatrix<double, GN::n_du_terms> du_at_gp;
    HybMatrix<double, GN::n_ddu_terms> ddu_at_gp;

    HybMatrix<double, GN::n_ddbath_terms> ddbath_at_gp;
    HybMatrix<double, GN::n_dddbath_terms> dddbath_at_gp;

    HybMatrix<double, GN::n_dimensions * GN::n_dimensions> w1_w1_kernel_at_gp;
    HybMatrix<double, GN::n_dimensions> w1_w2_kernel_at_gp;
    HybMatrix<double, GN::n_dimensions> w1_rhs_kernel_at_gp;

    HybMatrix<double, GN::n_dimensions> w2_w1_kernel_at_gp;
    DynRowVector<double> w2_w2_kernel_at_gp;

    DynMatrix<double> w1_w1;
    DynMatrix<double> w1_w2;
    DynVector<double> w1_rhs;

    DynMatrix<double> w2_w1;
    DynMatrix<double> w2_w2;
    DynMatrix<double> w2_w2_inv;
    /* rhs_w2 = 0 */
};
}
}

#endif
