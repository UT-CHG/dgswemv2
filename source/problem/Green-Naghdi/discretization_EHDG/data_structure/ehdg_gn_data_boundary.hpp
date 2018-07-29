#ifndef EHDG_GN_DATA_BOUNDARY_HPP
#define EHDG_GN_DATA_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : q_at_gp(GN::n_variables, ngp),
          aux_at_gp(GN::n_auxiliaries, ngp),
          Fn_at_gp(GN::n_variables, ngp),
          F_hat_at_gp(GN::n_variables, ngp) {}

    HybMatrix<double, GN::n_variables> q_at_gp;
    HybMatrix<double, GN::n_auxiliaries> aux_at_gp;

    HybMatrix<double, GN::n_variables> Fn_at_gp;
    HybMatrix<double, GN::n_variables> F_hat_at_gp;
};
}
}

#endif
