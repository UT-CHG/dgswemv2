#ifndef RKDG_SWE_PROC_BOUND_HPP
#define RKDG_SWE_PROC_BOUND_HPP

#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
template <typename BoundaryType>
void Problem::boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];
        auto& sp_at_gp = bound.data.spherical_projection.sp_at_gp_boundary[bound.bound_id];

        bound.ComputeUgp(state.ze, boundary.ze_at_gp);
        bound.ComputeUgp(state.qx, boundary.qx_at_gp);
        bound.ComputeUgp(state.qy, boundary.qy_at_gp);

        bound.boundary_condition.ComputeFlux(stepper,
                                             bound.surface_normal,
                                             sp_at_gp,
                                             boundary.bath_at_gp,
                                             boundary.ze_at_gp,
                                             boundary.qx_at_gp,
                                             boundary.qy_at_gp,
                                             boundary.ze_numerical_flux_at_gp,
                                             boundary.qx_numerical_flux_at_gp,
                                             boundary.qy_numerical_flux_at_gp);

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] -= bound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
            state.rhs_qx[dof] -= bound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
            state.rhs_qy[dof] -= bound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
        }
    }
}
}
}

#endif
