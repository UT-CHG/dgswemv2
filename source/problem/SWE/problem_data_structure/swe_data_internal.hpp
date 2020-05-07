#ifndef SWE_DATA_INTERNAL_HPP
#define SWE_DATA_INTERNAL_HPP

namespace SWE {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : q_at_gp(SWE::n_variables, ngp),
          aux_at_gp(SWE::n_auxiliaries, ngp),
          Fx_at_gp(SWE::n_variables, ngp),
          Fy_at_gp(SWE::n_variables, ngp),
          source_at_gp(SWE::n_variables, ngp),
          db_at_gp(SWE::n_dimensions, ngp),
          tau_s_at_gp(SWE::n_dimensions, ngp),
          dp_atm_at_gp(SWE::n_dimensions, ngp),
          dtide_pot_at_gp(SWE::n_dimensions, ngp),
          q_prev_at_gp(SWE::n_variables, ngp),
          del_q_DT_at_gp(SWE::n_variables, ngp),
          kronecker_DT_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          dFx_dq_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          dFy_dq_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          dsource_dq_at_gp(SWE::n_variables * SWE::n_variables, ngp) {}

    HybMatrix<double, SWE::n_variables> q_at_gp;
    HybMatrix<double, SWE::n_auxiliaries> aux_at_gp;

    HybMatrix<double, SWE::n_variables> Fx_at_gp;
    HybMatrix<double, SWE::n_variables> Fy_at_gp;

    HybMatrix<double, SWE::n_variables> source_at_gp;
    HybMatrix<double, SWE::n_dimensions> db_at_gp;
    HybMatrix<double, SWE::n_dimensions> tau_s_at_gp;
    HybMatrix<double, SWE::n_dimensions> dp_atm_at_gp;
    HybMatrix<double, SWE::n_dimensions> dtide_pot_at_gp;

    HybMatrix<double, SWE::n_variables> q_prev_at_gp;
    HybMatrix<double, SWE::n_variables> del_q_DT_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> kronecker_DT_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dFx_dq_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dFy_dq_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dsource_dq_at_gp;

    DynMatrix<double> delta_local_inv;
    DynMatrix<double> delta_local;
    DynVector<double> rhs_local;
    DynVector<double> rhs_prev;
};
}

#endif
