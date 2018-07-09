#ifndef EHDG_SWE_PROC_DBOUND_HPP
#define EHDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename DistributedBoundaryType>
void Problem::global_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.q, boundary.q_at_gp);

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        boundary.h_at_gp[gp] = boundary.q_at_gp[gp][SWE::Variables::ze] + boundary.bath_at_gp[gp];
    }

    /* Compute fluxes at boundary state */
    double nx, ny;
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        nx = dbound.surface_normal[gp][GlobalCoord::x];
        ny = dbound.surface_normal[gp][GlobalCoord::y];

        u = boundary.q_at_gp[gp][SWE::Variables::qx] / boundary.h_at_gp[gp];
        v = boundary.q_at_gp[gp][SWE::Variables::qy] / boundary.h_at_gp[gp];

        uuh = u * boundary.q_at_gp[gp][SWE::Variables::qx];
        vvh = v * boundary.q_at_gp[gp][SWE::Variables::qy];
        uvh = u * boundary.q_at_gp[gp][SWE::Variables::qy];
        pe  = Global::g * (0.5 * std::pow(boundary.q_at_gp[gp][SWE::Variables::ze], 2) +
                          boundary.q_at_gp[gp][SWE::Variables::ze] * boundary.bath_at_gp[gp]);

        boundary.Fn_at_gp[gp][SWE::Variables::ze] =
            boundary.q_at_gp[gp][SWE::Variables::qx] * nx + boundary.q_at_gp[gp][SWE::Variables::qy] * ny;
        boundary.Fn_at_gp[gp][SWE::Variables::qx] = (uuh + pe) * nx + uvh * ny;
        boundary.Fn_at_gp[gp][SWE::Variables::qy] = uvh * nx + (vvh + pe) * ny;
    }

    dbound.boundary_condition.exchanger.SetEX(boundary.q_at_gp, boundary.Fn_at_gp);
}

template <typename DistributedBoundaryType>
void Problem::local_distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < dbound.data.get_ndof(); ++dof) {
        state.rhs[dof] -= dbound.IntegrationPhi(dof, boundary.F_hat_at_gp);
    }
}
}
}

#endif
