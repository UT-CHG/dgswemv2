#ifndef IHDG_GN_PROC_BOUND_HPP
#define IHDG_GN_PROC_BOUND_HPP

namespace GN {
namespace IHDG {
template <typename BoundaryType>
void Problem::global_swe_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    boundary.q_at_gp = bound.ComputeUgp(state.q);

    row(boundary.aux_at_gp, GN::Auxiliaries::h) =
        row(boundary.q_at_gp, GN::Variables::ze) + row(boundary.aux_at_gp, GN::Auxiliaries::bath);

    /* Compute fluxes at boundary state */
    auto nx = row(bound.surface_normal, GlobalCoord::x);
    auto ny = row(bound.surface_normal, GlobalCoord::y);

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
}

template <typename BoundaryType>
void Problem::local_swe_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    // now compute contributions to the righthand side
    state.rhs -= bound.IntegrationPhi(boundary.F_hat_at_gp);
}
}
}

#endif
