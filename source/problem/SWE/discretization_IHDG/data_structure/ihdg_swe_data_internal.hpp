#ifndef IHDG_SWE_DATA_INTERNAL_HPP
#define IHDG_SWE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : q_at_gp(ngp),
          aux_at_gp(ngp),
          q_prev_at_gp(ngp),
          aux_prev_at_gp(ngp),
          Fx_at_gp(ngp),
          Fy_at_gp(ngp),
          source_at_gp(ngp),
          dbath_at_gp(ngp),
          tau_s_at_gp(ngp),
          dp_atm_at_gp(ngp),
          dtide_pot_at_gp(ngp),
          del_q_DT_at_gp(ngp),
          kronecker_DT_at_gp(ngp),
          dFx_dq_at_gp(ngp),
          dFy_dq_at_gp(ngp) {}

    std::vector<StatVector<double, SWE::n_variables>> q_at_gp;
    std::vector<StatVector<double, SWE::n_auxiliaries>> aux_at_gp;

    std::vector<StatVector<double, SWE::n_variables>> q_prev_at_gp;
    std::vector<StatVector<double, SWE::n_auxiliaries>> aux_prev_at_gp;

    std::vector<StatVector<double, SWE::n_variables>> Fx_at_gp;
    std::vector<StatVector<double, SWE::n_variables>> Fy_at_gp;

    std::vector<StatVector<double, SWE::n_variables>> source_at_gp;
    std::vector<StatVector<double, SWE::n_dimensions>> dbath_at_gp;
    std::vector<StatVector<double, SWE::n_dimensions>> tau_s_at_gp;
    std::vector<StatVector<double, SWE::n_dimensions>> dp_atm_at_gp;
    std::vector<StatVector<double, SWE::n_dimensions>> dtide_pot_at_gp;

    std::vector<StatVector<double, SWE::n_variables>> del_q_DT_at_gp;
    std::vector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> kronecker_DT_at_gp;
    std::vector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dFx_dq_at_gp;
    std::vector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dFy_dq_at_gp;

    DynMatrix<double> delta_local;
    DynVector<double> rhs_local;
};
}
}

#endif
