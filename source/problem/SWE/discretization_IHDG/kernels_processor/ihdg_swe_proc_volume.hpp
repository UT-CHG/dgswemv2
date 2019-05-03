#ifndef IHDG_SWE_PROC_VOLUME_HPP
#define IHDG_SWE_PROC_VOLUME_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace IHDG {
template <typename ElementType>
void Problem::init_volume_kernel(const ProblemStepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state_prev = elt.data.state[stage];
    auto& state_next = elt.data.state[stage + 1];
    auto& internal   = elt.data.internal;

    // Initial guess of state
    state_next.q = state_prev.q;

    internal.q_at_gp      = elt.ComputeUgp(state_prev.q);
    internal.q_prev_at_gp = internal.q_at_gp;

    set_constant(internal.rhs_prev, 0.0);

    if (stepper.GetOrder() == 2) {
        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        SWE::get_F(internal.q_at_gp, internal.aux_at_gp, internal.Fx_at_gp, internal.Fy_at_gp);

        for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
            subvector(internal.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) =
                elt.IntegrationDPhi(GlobalCoord::x, dof_i, internal.Fx_at_gp) +
                elt.IntegrationDPhi(GlobalCoord::y, dof_i, internal.Fy_at_gp);
        }
    }
}

template <typename ElementType>
void Problem::local_volume_kernel(const ProblemStepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage + 1];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    row(internal.aux_at_gp, SWE::Auxiliaries::h) =
        row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

    SWE::get_F(internal.q_at_gp, internal.aux_at_gp, internal.Fx_at_gp, internal.Fy_at_gp);

    // del_q / DT
    internal.del_q_DT_at_gp = (internal.q_at_gp - internal.q_prev_at_gp) / stepper.GetDT();

    SWE::get_dF_dq(internal.q_at_gp, internal.aux_at_gp, internal.dFx_dq_at_gp, internal.dFy_dq_at_gp);

    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        // Kronecker delta / DT
        column(internal.kronecker_DT_at_gp, gp) = IdentityVector<double>(SWE::n_variables) / stepper.GetDT();
    }

    double theta = stepper.GetTheta();

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    elt.IntegrationPhiPhi(dof_j, dof_i, internal.kronecker_DT_at_gp) -
                    (1.0 - theta) * elt.IntegrationPhiDPhi(dof_j, GlobalCoord::x, dof_i, internal.dFx_dq_at_gp) -
                    (1.0 - theta) * elt.IntegrationPhiDPhi(dof_j, GlobalCoord::y, dof_i, internal.dFy_dq_at_gp));
        }

        subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) =
            -elt.IntegrationPhi(dof_i, internal.del_q_DT_at_gp) +
            (1.0 - theta) * elt.IntegrationDPhi(GlobalCoord::x, dof_i, internal.Fx_at_gp) +
            (1.0 - theta) * elt.IntegrationDPhi(GlobalCoord::y, dof_i, internal.Fy_at_gp);
    }

    internal.rhs_local += theta * internal.rhs_prev;
}
}
}

#endif
