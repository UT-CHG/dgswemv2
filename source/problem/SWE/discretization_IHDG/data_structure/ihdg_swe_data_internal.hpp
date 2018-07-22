#ifndef IHDG_SWE_DATA_INTERNAL_HPP
#define IHDG_SWE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
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
          dtide_pot_at_gp(SWE::n_dimensions, ngp),
          q_prev_at_gp(SWE::n_variables, ngp),
          del_q_DT_at_gp(SWE::n_variables, ngp),
          kronecker_DT_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          dFx_dq_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          dFy_dq_at_gp(SWE::n_variables * SWE::n_variables, ngp) {}

    DynMatrix<double> q_at_gp;
    DynMatrix<double> aux_at_gp;

    DynMatrix<double> Fx_at_gp;
    DynMatrix<double> Fy_at_gp;

    DynMatrix<double> source_at_gp;
    DynMatrix<double> dbath_at_gp;
    DynMatrix<double> tau_s_at_gp;
    DynMatrix<double> dp_atm_at_gp;
    DynMatrix<double> dtide_pot_at_gp;

    DynMatrix<double> q_prev_at_gp;
    DynMatrix<double> del_q_DT_at_gp;
    DynMatrix<double> kronecker_DT_at_gp;
    DynMatrix<double> dFx_dq_at_gp;
    DynMatrix<double> dFy_dq_at_gp;

    DynMatrix<double> delta_local;
    DynVector<double> rhs_local;
};
}
}

#endif
