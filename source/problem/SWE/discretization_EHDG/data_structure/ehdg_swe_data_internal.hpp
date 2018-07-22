#ifndef EHDG_SWE_DATA_INTERNAL_HPP
#define EHDG_SWE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : q_at_gp(SWE::n_variables, ngp),
          aux_at_gp(SWE::n_auxiliaries, ngp),
          Fx_at_gp(SWE::n_variables, ngp),
          Fy_at_gp(SWE::n_variables, ngp),
          source_at_gp(SWE::n_variables, ngp),
          dbath_at_gp(SWE::n_dimensions, ngp),
          tau_s_at_gp(SWE::n_dimensions, ngp),
          dp_atm_at_gp(SWE::n_dimensions, ngp),
          dtide_pot_at_gp(SWE::n_dimensions, ngp) {}

    DynMatrix<double> q_at_gp;
    DynMatrix<double> aux_at_gp;

    DynMatrix<double> Fx_at_gp;
    DynMatrix<double> Fy_at_gp;

    DynMatrix<double> source_at_gp;
    DynMatrix<double> dbath_at_gp;
    DynMatrix<double> tau_s_at_gp;
    DynMatrix<double> dp_atm_at_gp;
    DynMatrix<double> dtide_pot_at_gp;
};
}
}

#endif
