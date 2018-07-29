#ifndef EHDG_GN_DATA_DISP_CORR_HPP
#define EHDG_GN_DATA_DISP_CORR_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct DispersiveCorrection {
    DispersiveCorrection() = default;
    DispersiveCorrection(const uint ndof, const uint ngp_internal, const std::vector<uint>& ngp_boundary)
        : du(GN::n_du_terms, ndof),
          ddu(GN::n_ddu_terms, ndof),
          u_at_gp(GN::n_dimensions, ngp_internal),
          du_at_gp(GN::n_du_terms, ngp_internal),
          u_hat_at_gp(ngp_boundary.size()),
          du_hat_at_gp(ngp_boundary.size()) {
        // *** //
        for (uint bound = 0; bound < ngp_boundary.size(); ++bound) {
            this->u_hat_at_gp[bound].resize(GN::n_dimensions, ngp_boundary[bound]);
            this->du_hat_at_gp[bound].resize(GN::n_du_terms, ngp_boundary[bound]);
        }
    }

    HybMatrix<double, GN::n_du_terms> du;
    HybMatrix<double, GN::n_ddu_terms> ddu;

    HybMatrix<double, GN::n_dimensions> u_at_gp;
    HybMatrix<double, GN::n_du_terms> du_at_gp;

    std::vector<HybMatrix<double, GN::n_dimensions>> u_hat_at_gp;
    std::vector<HybMatrix<double, GN::n_du_terms>> du_hat_at_gp;
};
}
}

#endif