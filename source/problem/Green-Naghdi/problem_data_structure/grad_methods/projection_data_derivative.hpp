#ifndef PROJECTION_DATA_DERIVATIVE_HPP
#define PROJECTION_DATA_DERIVATIVE_HPP

namespace GN {
struct Derivative {
    Derivative() = default;
    Derivative(const uint nvrtx, const uint nbound, const std::vector<uint>& ngp_boundary)
        : bath_hat_at_gp(nbound),
          ddbath_hat_at_gp(nbound),
          ze_hat_at_gp(nbound),
          u_hat_at_gp(nbound),
          du_hat_at_gp(nbound) {
        for (uint bound = 0; bound < nbound; ++bound) {
            bath_hat_at_gp[bound]   = DynRowVector<double>(ngp_boundary[bound]);
            ddbath_hat_at_gp[bound] = HybMatrix<double, GN::n_ddbath_terms>(GN::n_ddbath_terms, ngp_boundary[bound]);

            bath_hat_at_gp[bound] = DynRowVector<double>(ngp_boundary[bound]);
            u_hat_at_gp[bound]    = HybMatrix<double, GN::n_dimensions>(GN::n_dimensions, ngp_boundary[bound]);
            du_hat_at_gp[bound]   = HybMatrix<double, GN::n_du_terms>(GN::n_du_terms, ngp_boundary[bound]);
        }
    }

    std::vector<DynRowVector<double>> bath_hat_at_gp;
    AlignedVector<HybMatrix<double, GN::n_ddbath_terms>> ddbath_hat_at_gp;

    std::vector<DynRowVector<double>> ze_hat_at_gp;
    AlignedVector<HybMatrix<double, GN::n_dimensions>> u_hat_at_gp;
    AlignedVector<HybMatrix<double, GN::n_du_terms>> du_hat_at_gp;
};
}

#endif
