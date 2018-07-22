#ifndef EHDG_SWE_DATA_BOUNDARY_HPP
#define EHDG_SWE_DATA_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : q_at_gp(SWE::n_variables, ngp),
          aux_at_gp(SWE::n_auxiliaries, ngp),
          Fn_at_gp(SWE::n_variables, ngp),
          F_hat_at_gp(SWE::n_variables, ngp) {}

    DynMatrix<double> q_at_gp;
    DynMatrix<double> aux_at_gp;

    DynMatrix<double> Fn_at_gp;
    DynMatrix<double> F_hat_at_gp;
};
}
}

#endif
