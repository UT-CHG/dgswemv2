#ifndef EHDG_SWE_DATA_BOUNDARY_HPP
#define EHDG_SWE_DATA_BOUNDARY_HPP

namespace SWE {
namespace EHDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : q_at_gp(SWE::n_variables, ngp),
          aux_at_gp(SWE::n_auxiliaries, ngp),
          Fn_at_gp(SWE::n_variables, ngp),
          F_hat_at_gp(SWE::n_variables, ngp) {}

    HybMatrix<double, SWE::n_variables> q_at_gp;
    HybMatrix<double, SWE::n_auxiliaries> aux_at_gp;

    HybMatrix<double, SWE::n_variables> Fn_at_gp;
    HybMatrix<double, SWE::n_variables> F_hat_at_gp;
};
}
}

#endif
