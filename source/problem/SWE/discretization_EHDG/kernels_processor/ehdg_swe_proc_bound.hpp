#ifndef EHDG_SWE_PROC_BOUND_HPP
#define EHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace EHDG {
template <typename BoundaryType>
void Problem::global_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    boundary.q_at_gp = bound.ComputeUgp(state.q);

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        boundary.aux_at_gp(SWE::Auxiliaries::h, gp) =
            boundary.q_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
    }

    /* Compute fluxes at boundary state */
    double nx, ny;
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        nx = bound.surface_normal(GlobalCoord::x, gp);
        ny = bound.surface_normal(GlobalCoord::y, gp);

        u = boundary.q_at_gp(SWE::Variables::qx, gp) / boundary.aux_at_gp(SWE::Auxiliaries::h, gp);
        v = boundary.q_at_gp(SWE::Variables::qy, gp) / boundary.aux_at_gp(SWE::Auxiliaries::h, gp);

        uuh = u * boundary.q_at_gp(SWE::Variables::qx, gp);
        vvh = v * boundary.q_at_gp(SWE::Variables::qy, gp);
        uvh = u * boundary.q_at_gp(SWE::Variables::qy, gp);
        pe  = Global::g * (0.5 * std::pow(boundary.q_at_gp(SWE::Variables::ze, gp), 2) +
                          boundary.q_at_gp(SWE::Variables::ze, gp) * boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));

        boundary.Fn_at_gp(SWE::Variables::ze, gp) =
            boundary.q_at_gp(SWE::Variables::qx, gp) * nx + boundary.q_at_gp(SWE::Variables::qy, gp) * ny;
        boundary.Fn_at_gp(SWE::Variables::qx, gp) = (uuh + pe) * nx + uvh * ny;
        boundary.Fn_at_gp(SWE::Variables::qy, gp) = uvh * nx + (vvh + pe) * ny;
    }
}

template <typename BoundaryType>
void Problem::local_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
        column(state.rhs, dof) -= bound.IntegrationPhi(dof, boundary.F_hat_at_gp);
    }
}
}
}

#endif
