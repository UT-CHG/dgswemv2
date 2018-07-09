#ifndef IHDG_SWE_PROC_BOUND_HPP
#define IHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace IHDG {
template <typename BoundaryType>
void Problem::prepare_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    bound.ComputeUgp(state.ze, boundary.ze_at_gp);
    bound.ComputeUgp(state.qx, boundary.qx_at_gp);
    bound.ComputeUgp(state.qy, boundary.qy_at_gp);

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        boundary.aux_at_gp[gp][SWE::Auxiliaries::h] =
            boundary.ze_at_gp[gp] + boundary.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    /* Compute fluxes at boundary state */
    double nx, ny;
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        nx = bound.surface_normal[gp][GlobalCoord::x];
        ny = bound.surface_normal[gp][GlobalCoord::y];

        u = boundary.qx_at_gp[gp] / boundary.aux_at_gp[gp][SWE::Auxiliaries::h];
        v = boundary.qy_at_gp[gp] / boundary.aux_at_gp[gp][SWE::Auxiliaries::h];

        uuh = u * boundary.qx_at_gp[gp];
        vvh = v * boundary.qy_at_gp[gp];
        uvh = u * boundary.qy_at_gp[gp];
        pe  = Global::g * (0.5 * std::pow(boundary.ze_at_gp[gp], 2) +
                          boundary.ze_at_gp[gp] * boundary.aux_at_gp[gp][SWE::Auxiliaries::bath]);

        boundary.ze_flux_dot_n_at_gp[gp] = boundary.qx_at_gp[gp] * nx + boundary.qy_at_gp[gp] * ny;
        boundary.qx_flux_dot_n_at_gp[gp] = (uuh + pe) * nx + uvh * ny;
        boundary.qy_flux_dot_n_at_gp[gp] = uvh * nx + (vvh + pe) * ny;
    }
}
}
}

#endif
