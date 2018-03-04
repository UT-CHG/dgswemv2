#ifndef SWE_PROC_BOUND_HPP
#define SWE_PROC_BOUND_HPP

#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
template <typename BoundaryType>
void Problem::boundary_kernel(const Stepper& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.get_stage();

        auto& state = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        bound.ComputeUgp(state.ze, boundary.ze_at_gp);
        bound.ComputeUgp(state.qx, boundary.qx_at_gp);
        bound.ComputeUgp(state.qy, boundary.qy_at_gp);

        double ze_ex, qx_ex, qy_ex;
        for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
            bound.boundary_condition.GetEX(stepper, gp, bound.surface_normal, boundary.ze_at_gp, boundary.qx_at_gp,
                                           boundary.qy_at_gp, ze_ex, qx_ex, qy_ex);

            LLF_flux(boundary.ze_at_gp[gp], ze_ex, boundary.qx_at_gp[gp], qx_ex, boundary.qy_at_gp[gp], qy_ex,
                     boundary.bath_at_gp[gp], bound.surface_normal[gp], boundary.ze_numerical_flux_at_gp[gp],
                     boundary.qx_numerical_flux_at_gp[gp], boundary.qy_numerical_flux_at_gp[gp]);
        }

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
