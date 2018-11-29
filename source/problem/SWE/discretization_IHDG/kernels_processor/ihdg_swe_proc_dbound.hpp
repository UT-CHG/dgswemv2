#ifndef IHDG_SWE_PROC_DBOUND_HPP
#define IHDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename DistributedBoundaryType>
void Problem::init_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound) {
    if (stepper.GetOrder() == 2) {
        const uint stage = stepper.GetStage();

        auto& state_prev = dbound.data.state[stage];
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        boundary.q_at_gp = dbound.ComputeUgp(state_prev.q);

        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        /* Compute fluxes at boundary state */
        auto nx = row(dbound.surface_normal, GlobalCoord::x);
        auto ny = row(dbound.surface_normal, GlobalCoord::y);

        auto u = vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
        auto v = vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

        auto uuh = vec_cw_mult(u, row(boundary.q_at_gp, SWE::Variables::qx));
        auto vvh = vec_cw_mult(v, row(boundary.q_at_gp, SWE::Variables::qy));
        auto uvh = vec_cw_mult(u, row(boundary.q_at_gp, SWE::Variables::qy));
        auto pe =
            Global::g *
            (0.5 * vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(boundary.q_at_gp, SWE::Variables::ze)) +
             vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

        // Fn terms
        row(boundary.F_hat_at_gp, SWE::Variables::ze) = vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qx), nx) +
                                                        vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qy), ny);
        row(boundary.F_hat_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
        row(boundary.F_hat_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);
    }
}

template <typename StepperType, typename DistributedBoundaryType>
void Problem::local_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage + 1];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    boundary.q_at_gp = dbound.ComputeUgp(state.q);

    row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
        row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary state */
    auto nx = row(dbound.surface_normal, GlobalCoord::x);
    auto ny = row(dbound.surface_normal, GlobalCoord::y);

    auto u = vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(boundary.q_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(boundary.q_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(boundary.q_at_gp, SWE::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(boundary.q_at_gp, SWE::Variables::ze)) +
               vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    // Fn terms
    row(boundary.F_hat_at_gp, SWE::Variables::ze) = vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qx), nx) +
                                                    vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qy), ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

    // dFn/dq terms
    set_constant(row(boundary.dF_hat_dq_at_gp, JacobianVariables::ze_ze), 0.0);
    row(boundary.dF_hat_dq_at_gp, JacobianVariables::ze_qx) = nx;
    row(boundary.dF_hat_dq_at_gp, JacobianVariables::ze_qy) = ny;

    row(boundary.dF_hat_dq_at_gp, JacobianVariables::qx_ze) =
        vec_cw_mult(-vec_cw_mult(u, u) + Global::g * row(boundary.aux_at_gp, SWE::Auxiliaries::h), nx) -
        vec_cw_mult(vec_cw_mult(u, v), ny);
    row(boundary.dF_hat_dq_at_gp, JacobianVariables::qx_qx) = 2.0 * vec_cw_mult(u, nx) + vec_cw_mult(v, ny);
    row(boundary.dF_hat_dq_at_gp, JacobianVariables::qx_qy) = vec_cw_mult(u, ny);

    row(boundary.dF_hat_dq_at_gp, JacobianVariables::qy_ze) =
        -vec_cw_mult(vec_cw_mult(u, v), nx) +
        vec_cw_mult(-vec_cw_mult(v, v) + Global::g * row(boundary.aux_at_gp, SWE::Auxiliaries::h), ny);
    row(boundary.dF_hat_dq_at_gp, JacobianVariables::qy_qx) = vec_cw_mult(v, nx);
    row(boundary.dF_hat_dq_at_gp, JacobianVariables::qy_qy) = vec_cw_mult(u, nx) + 2.0 * vec_cw_mult(v, ny);
}
}
}

#endif
