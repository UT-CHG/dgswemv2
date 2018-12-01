#ifndef IHDG_SWE_PROC_VOLUME_HPP
#define IHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename ElementType>
void Problem::init_volume_kernel(const StepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state_prev = elt.data.state[stage];
    auto& state_next = elt.data.state[stage + 1];
    auto& internal   = elt.data.internal;

    // Initial guess of state
    state_next.q = state_prev.q;

    internal.q_at_gp      = elt.ComputeUgp(state_prev.q);
    internal.q_prev_at_gp = internal.q_at_gp;

    if (stepper.GetOrder() == 2) {
        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        auto u = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
        auto v = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

        auto uuh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qx));
        auto vvh = vec_cw_mult(v, row(internal.q_at_gp, SWE::Variables::qy));
        auto uvh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qy));
        auto pe =
            Global::g *
            (0.5 * vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.q_at_gp, SWE::Variables::ze)) +
             vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.aux_at_gp, SWE::Auxiliaries::bath)));

        // Flux terms
        row(internal.Fx_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qx);
        row(internal.Fx_at_gp, SWE::Variables::qx) = uuh + pe;
        row(internal.Fx_at_gp, SWE::Variables::qy) = uvh;

        row(internal.Fy_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qy);
        row(internal.Fy_at_gp, SWE::Variables::qx) = uvh;
        row(internal.Fy_at_gp, SWE::Variables::qy) = vvh + pe;

        double theta = stepper.GetTheta();

        for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
            subvector(internal.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) =
                theta * elt.IntegrationDPhi(GlobalCoord::x, dof_i, internal.Fx_at_gp) +
                theta * elt.IntegrationDPhi(GlobalCoord::y, dof_i, internal.Fy_at_gp);
        }
    }
}

template <typename StepperType, typename ElementType>
void Problem::local_volume_kernel(const StepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage + 1];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    row(internal.aux_at_gp, SWE::Auxiliaries::h) =
        row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

    auto u = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(internal.q_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.q_at_gp, SWE::Variables::ze)) +
               vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.aux_at_gp, SWE::Auxiliaries::bath)));

    // Flux terms
    row(internal.Fx_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qx);
    row(internal.Fx_at_gp, SWE::Variables::qx) = uuh + pe;
    row(internal.Fx_at_gp, SWE::Variables::qy) = uvh;

    row(internal.Fy_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qy);
    row(internal.Fy_at_gp, SWE::Variables::qx) = uvh;
    row(internal.Fy_at_gp, SWE::Variables::qy) = vvh + pe;

    // del_q / DT
    internal.del_q_DT_at_gp = (internal.q_at_gp - internal.q_prev_at_gp) / stepper.GetDT();

    // dFx/dq terms
    set_constant(row(internal.dFx_dq_at_gp, JacobianVariables::ze_ze), 0.0);
    set_constant(row(internal.dFx_dq_at_gp, JacobianVariables::ze_qx), 1.0);
    set_constant(row(internal.dFx_dq_at_gp, JacobianVariables::ze_qy), 0.0);

    row(internal.dFx_dq_at_gp, JacobianVariables::qx_ze) =
        -vec_cw_mult(u, u) + Global::g * row(internal.aux_at_gp, SWE::Auxiliaries::h);
    row(internal.dFx_dq_at_gp, JacobianVariables::qx_qx) = 2.0 * u;
    set_constant(row(internal.dFx_dq_at_gp, JacobianVariables::qx_qy), 0.0);

    row(internal.dFx_dq_at_gp, JacobianVariables::qy_ze) = -vec_cw_mult(u, v);
    row(internal.dFx_dq_at_gp, JacobianVariables::qy_qx) = v;
    row(internal.dFx_dq_at_gp, JacobianVariables::qy_qy) = u;

    // dFy/dq terms
    set_constant(row(internal.dFy_dq_at_gp, JacobianVariables::ze_ze), 0.0);
    set_constant(row(internal.dFy_dq_at_gp, JacobianVariables::ze_qx), 0.0);
    set_constant(row(internal.dFy_dq_at_gp, JacobianVariables::ze_qy), 1.0);

    row(internal.dFy_dq_at_gp, JacobianVariables::qx_ze) = -vec_cw_mult(u, v);
    row(internal.dFy_dq_at_gp, JacobianVariables::qx_qx) = v;
    row(internal.dFy_dq_at_gp, JacobianVariables::qx_qy) = u;

    row(internal.dFy_dq_at_gp, JacobianVariables::qy_ze) =
        -vec_cw_mult(v, v) + Global::g * row(internal.aux_at_gp, SWE::Auxiliaries::h);
    set_constant(row(internal.dFy_dq_at_gp, JacobianVariables::qy_qx), 0.0);
    row(internal.dFy_dq_at_gp, JacobianVariables::qy_qy) = 2.0 * v;

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
