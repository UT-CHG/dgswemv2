#ifndef EHDG_GN_PROC_DBOUND_HPP
#define EHDG_GN_PROC_DBOUND_HPP

namespace GN {
namespace EHDG {
template <typename DistributedBoundaryType>
void Problem::global_swe_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    boundary.q_at_gp = dbound.ComputeUgp(state.q);

    row(boundary.aux_at_gp, GN::Auxiliaries::h) =
        row(boundary.q_at_gp, GN::Variables::ze) + row(boundary.aux_at_gp, GN::Auxiliaries::bath);

    /* Compute fluxes at boundary state */
    auto nx = row(dbound.surface_normal, GlobalCoord::x);
    auto ny = row(dbound.surface_normal, GlobalCoord::y);

    auto u = vec_cw_div(row(boundary.q_at_gp, GN::Variables::qx), row(boundary.aux_at_gp, GN::Auxiliaries::h));
    auto v = vec_cw_div(row(boundary.q_at_gp, GN::Variables::qy), row(boundary.aux_at_gp, GN::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(boundary.q_at_gp, GN::Variables::qx));
    auto vvh = vec_cw_mult(v, row(boundary.q_at_gp, GN::Variables::qy));
    auto uvh = vec_cw_mult(u, row(boundary.q_at_gp, GN::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(boundary.q_at_gp, GN::Variables::ze), row(boundary.q_at_gp, GN::Variables::ze)) +
               vec_cw_mult(row(boundary.q_at_gp, GN::Variables::ze), row(boundary.aux_at_gp, GN::Auxiliaries::bath)));

    row(boundary.Fn_at_gp, GN::Variables::ze) = vec_cw_mult(row(boundary.q_at_gp, GN::Variables::qx), nx) +
                                                vec_cw_mult(row(boundary.q_at_gp, GN::Variables::qy), ny);
    row(boundary.Fn_at_gp, GN::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
    row(boundary.Fn_at_gp, GN::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

    dbound.boundary_condition.exchanger.SetEX(boundary.q_at_gp, boundary.Fn_at_gp);
}

template <typename DistributedBoundaryType>
void Problem::local_swe_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    // now compute contributions to the righthand side
    state.rhs -= dbound.IntegrationPhi(boundary.F_hat_at_gp);
}
}
}

#endif
