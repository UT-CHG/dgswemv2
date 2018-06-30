#ifndef EHDG_SWE_PROC_BOUND_HPP
#define EHDG_SWE_PROC_BOUND_HPP

#include "problem/SWE/discretization_EHDG/numerical_fluxes/ehdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace EHDG {
template <typename BoundaryType>
void Problem::global_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    bound.ComputeUgp(state.ze, boundary.ze_at_gp);
    bound.ComputeUgp(state.qx, boundary.qx_at_gp);
    bound.ComputeUgp(state.qy, boundary.qy_at_gp);

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        boundary.h_at_gp[gp] = boundary.ze_at_gp[gp] + boundary.bath_at_gp[gp];
    }

    /* Compute fluxes at boundary state */
    double u, v;
    double uuh, vvh, uvh, pe;

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        u = boundary.qx_at_gp[gp] / boundary.h_at_gp[gp];
        v = boundary.qy_at_gp[gp] / boundary.h_at_gp[gp];

        uuh = u * boundary.qx_at_gp[gp];
        vvh = v * boundary.qy_at_gp[gp];
        uvh = u * boundary.qy_at_gp[gp];
        pe  = Global::g * (0.5 * std::pow(boundary.ze_at_gp[gp], 2) + boundary.ze_at_gp[gp] * boundary.bath_at_gp[gp]);

        boundary.ze_flux_at_gp[GlobalCoord::x][gp] = boundary.qx_at_gp[gp];
        boundary.ze_flux_at_gp[GlobalCoord::y][gp] = boundary.qy_at_gp[gp];

        boundary.qx_flux_at_gp[GlobalCoord::x][gp] = (uuh + pe);
        boundary.qx_flux_at_gp[GlobalCoord::y][gp] = uvh;

        boundary.qy_flux_at_gp[GlobalCoord::x][gp] = uvh;
        boundary.qy_flux_at_gp[GlobalCoord::y][gp] = vvh + pe;
    }
}

template <typename BoundaryType>
void Problem::local_boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] -= bound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
        state.rhs_qx[dof] -= bound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
        state.rhs_qy[dof] -= bound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
    }
}
}
}

#endif
