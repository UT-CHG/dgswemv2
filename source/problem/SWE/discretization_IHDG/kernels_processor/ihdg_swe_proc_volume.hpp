#ifndef IHDG_SWE_PROC_VOLUME_HPP
#define IHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace IHDG {
template <typename ElementType>
void Problem::local_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state_prev = elt.data.state[stage];
    auto& state      = elt.data.state[stage + 1];
    auto& internal   = elt.data.internal;

    elt.ComputeUgp(state.q, internal.q_at_gp);
    elt.ComputeUgp(state_prev.q, internal.q_prev_at_gp);

    double u = 0.0;
    double v = 0.0;

    double uuh = 0.0;
    double vvh = 0.0;
    double uvh = 0.0;
    double pe  = 0.0;

    // assemble flux
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.aux_at_gp[gp][SWE::Auxiliaries::h] =
            internal.q_at_gp[gp][SWE::Variables::ze] + internal.aux_at_gp[gp][SWE::Auxiliaries::bath];

        u = internal.q_at_gp[gp][SWE::Variables::qx] / internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        v = internal.q_at_gp[gp][SWE::Variables::qy] / internal.aux_at_gp[gp][SWE::Auxiliaries::h];

        uuh = u * internal.q_at_gp[gp][SWE::Variables::qx];
        vvh = v * internal.q_at_gp[gp][SWE::Variables::qy];
        uvh = u * internal.q_at_gp[gp][SWE::Variables::qy];
        pe  = Global::g * (0.5 * std::pow(internal.q_at_gp[gp][SWE::Variables::ze], 2) +
                          internal.q_at_gp[gp][SWE::Variables::ze] * internal.aux_at_gp[gp][SWE::Auxiliaries::bath]);

        // Flux terms
        internal.Fx_at_gp[gp][SWE::Variables::ze] = internal.q_at_gp[gp][SWE::Variables::qx];
        internal.Fx_at_gp[gp][SWE::Variables::qx] = uuh + pe;
        internal.Fx_at_gp[gp][SWE::Variables::qy] = uvh;

        internal.Fy_at_gp[gp][SWE::Variables::ze] = internal.q_at_gp[gp][SWE::Variables::qy];
        internal.Fy_at_gp[gp][SWE::Variables::qx] = uvh;
        internal.Fy_at_gp[gp][SWE::Variables::qy] = vvh + pe;

        // del_q / DT
        internal.del_q_DT_at_gp[gp] = (internal.q_at_gp[gp] - internal.q_prev_at_gp[gp]) / stepper.GetDT();

        // Kronecker delta / DT
        internal.kronecker_DT_at_gp[gp] = blaze::IdentityMatrix<double>(SWE::n_variables) / stepper.GetDT();

        // dFx/dq terms
        internal.dFx_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        internal.dFx_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qx) = 1.0;
        internal.dFx_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qy) = 0.0;

        internal.dFx_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::ze) =
            -u * u + Global::g * internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        internal.dFx_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qx) = 2 * u;
        internal.dFx_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qy) = 0.0;

        internal.dFx_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::ze) = -u * v;
        internal.dFx_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qx) = v;
        internal.dFx_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qy) = u;

        // dFy/dq terms
        internal.dFy_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        internal.dFy_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qx) = 0.0;
        internal.dFy_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qy) = 1.0;

        internal.dFy_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::ze) = -u * v;
        internal.dFy_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qx) = v;
        internal.dFy_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qy) = u;

        internal.dFy_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::ze) =
            -v * v + Global::g * internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        internal.dFy_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qx) = 0.0;
        internal.dFy_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qy) = 2 * v;
    }

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); dof_j++) {
            blaze::submatrix(internal.delta_local,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) =
                elt.IntegrationPhiPhi(dof_i, dof_j, internal.kronecker_DT_at_gp) -
                elt.IntegrationPhiDPhi(dof_j, GlobalCoord::x, dof_i, internal.dFx_dq_at_gp) -
                elt.IntegrationPhiDPhi(dof_j, GlobalCoord::y, dof_i, internal.dFy_dq_at_gp);
        }

        blaze::subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) =
            -elt.IntegrationPhi(dof_i, internal.del_q_DT_at_gp) +
            elt.IntegrationDPhi(GlobalCoord::x, dof_i, internal.Fx_at_gp) +
            elt.IntegrationDPhi(GlobalCoord::y, dof_i, internal.Fy_at_gp);
    }
}
}
}

#endif
