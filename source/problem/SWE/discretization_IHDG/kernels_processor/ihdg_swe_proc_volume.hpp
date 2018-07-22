#ifndef IHDG_SWE_PROC_VOLUME_HPP
#define IHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace IHDG {
template <typename ElementType>
void Problem::local_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage + 1];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    double u = 0.0;
    double v = 0.0;

    double uuh = 0.0;
    double vvh = 0.0;
    double uvh = 0.0;
    double pe  = 0.0;

    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.aux_at_gp(SWE::Auxiliaries::h, gp) =
            internal.q_at_gp(SWE::Variables::ze, gp) + internal.aux_at_gp(SWE::Auxiliaries::bath, gp);

        u = internal.q_at_gp(SWE::Variables::qx, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);
        v = internal.q_at_gp(SWE::Variables::qy, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);

        uuh = u * internal.q_at_gp(SWE::Variables::qx, gp);
        vvh = v * internal.q_at_gp(SWE::Variables::qy, gp);
        uvh = u * internal.q_at_gp(SWE::Variables::qy, gp);
        pe  = Global::g * (0.5 * std::pow(internal.q_at_gp(SWE::Variables::ze, gp), 2) +
                          internal.q_at_gp(SWE::Variables::ze, gp) * internal.aux_at_gp(SWE::Auxiliaries::bath, gp));

        // Flux terms
        internal.Fx_at_gp(SWE::Variables::ze, gp) = internal.q_at_gp(SWE::Variables::qx, gp);
        internal.Fx_at_gp(SWE::Variables::qx, gp) = uuh + pe;
        internal.Fx_at_gp(SWE::Variables::qy, gp) = uvh;

        internal.Fy_at_gp(SWE::Variables::ze, gp) = internal.q_at_gp(SWE::Variables::qy, gp);
        internal.Fy_at_gp(SWE::Variables::qx, gp) = uvh;
        internal.Fy_at_gp(SWE::Variables::qy, gp) = vvh + pe;

        // del_q / DT
        column(internal.del_q_DT_at_gp, gp) =
            (column(internal.q_at_gp, gp) - column(internal.q_prev_at_gp, gp)) / stepper.GetDT();

        // Kronecker delta / DT
        column(internal.kronecker_DT_at_gp, gp) = I_vector / stepper.GetDT();

        // dFx/dq terms
        internal.dFx_dq_at_gp(JacobianVariables::ze_ze, gp) = 0.0;
        internal.dFx_dq_at_gp(JacobianVariables::ze_qx, gp) = 1.0;
        internal.dFx_dq_at_gp(JacobianVariables::ze_qy, gp) = 0.0;

        internal.dFx_dq_at_gp(JacobianVariables::qx_ze, gp) =
            -u * u + Global::g * internal.aux_at_gp(SWE::Auxiliaries::h, gp);
        internal.dFx_dq_at_gp(JacobianVariables::qx_qx, gp) = 2 * u;
        internal.dFx_dq_at_gp(JacobianVariables::qx_qy, gp) = 0.0;

        internal.dFx_dq_at_gp(JacobianVariables::qy_ze, gp) = -u * v;
        internal.dFx_dq_at_gp(JacobianVariables::qy_qx, gp) = v;
        internal.dFx_dq_at_gp(JacobianVariables::qy_qy, gp) = u;

        // dFy/dq terms
        internal.dFy_dq_at_gp(JacobianVariables::ze_ze, gp) = 0.0;
        internal.dFy_dq_at_gp(JacobianVariables::ze_qx, gp) = 0.0;
        internal.dFy_dq_at_gp(JacobianVariables::ze_qy, gp) = 1.0;

        internal.dFy_dq_at_gp(JacobianVariables::qx_ze, gp) = -u * v;
        internal.dFy_dq_at_gp(JacobianVariables::qx_qx, gp) = v;
        internal.dFy_dq_at_gp(JacobianVariables::qx_qy, gp) = u;

        internal.dFy_dq_at_gp(JacobianVariables::qy_ze, gp) =
            -v * v + Global::g * internal.aux_at_gp(SWE::Auxiliaries::h, gp);
        internal.dFy_dq_at_gp(JacobianVariables::qy_qx, gp) = 0.0;
        internal.dFy_dq_at_gp(JacobianVariables::qy_qy, gp) = 2 * v;
    }

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); dof_j++) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    elt.IntegrationPhiPhi(dof_j, dof_i, internal.kronecker_DT_at_gp) -
                    elt.IntegrationPhiDPhi(dof_j, GlobalCoord::x, dof_i, internal.dFx_dq_at_gp) -
                    elt.IntegrationPhiDPhi(dof_j, GlobalCoord::y, dof_i, internal.dFy_dq_at_gp));
        }

        subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) =
            -elt.IntegrationPhi(dof_i, internal.del_q_DT_at_gp) +
            elt.IntegrationDPhi(GlobalCoord::x, dof_i, internal.Fx_at_gp) +
            elt.IntegrationDPhi(GlobalCoord::y, dof_i, internal.Fy_at_gp);
    }
}
}
}

#endif
