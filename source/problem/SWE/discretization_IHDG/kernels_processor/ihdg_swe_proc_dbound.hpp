#ifndef IHDG_SWE_PROC_DBOUND_HPP
#define IHDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace IHDG {
template <typename DistributedBoundaryType>
void Problem::local_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.ze, boundary.ze_at_gp);
    dbound.ComputeUgp(state.qx, boundary.qx_at_gp);
    dbound.ComputeUgp(state.qy, boundary.qy_at_gp);

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        boundary.aux_at_gp(SWE::Auxiliaries::h, gp) =
            boundary.ze_at_gp[gp] + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
    }

    /* Compute fluxes at boundary state */
    double nx, ny;
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        nx = dbound.surface_normal[gp][GlobalCoord::x];
        ny = dbound.surface_normal[gp][GlobalCoord::y];

        u = boundary.qx_at_gp[gp] / boundary.aux_at_gp(SWE::Auxiliaries::h, gp);
        v = boundary.qy_at_gp[gp] / boundary.aux_at_gp(SWE::Auxiliaries::h, gp);

        uuh = u * boundary.qx_at_gp[gp];
        vvh = v * boundary.qy_at_gp[gp];
        uvh = u * boundary.qy_at_gp[gp];
        pe  = Global::g * (0.5 * std::pow(boundary.ze_at_gp[gp], 2) +
                          boundary.ze_at_gp[gp] * boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));

        boundary.ze_flux_dot_n_at_gp[gp] = boundary.qx_at_gp[gp] * nx + boundary.qy_at_gp[gp] * ny;
        boundary.qx_flux_dot_n_at_gp[gp] = (uuh + pe) * nx + uvh * ny;
        boundary.qy_flux_dot_n_at_gp[gp] = uvh * nx + (vvh + pe) * ny;
    }

    dbound.boundary_condition.exchanger.SetEX(boundary.ze_at_gp,
                                              boundary.qx_at_gp,
                                              boundary.qy_at_gp,
                                              boundary.ze_flux_dot_n_at_gp,
                                              boundary.qx_flux_dot_n_at_gp,
                                              boundary.qy_flux_dot_n_at_gp);
}
}
}

#endif
