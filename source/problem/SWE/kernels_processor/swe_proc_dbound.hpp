#ifndef SWE_PROC_DBOUND_HPP
#define SWE_PROC_DBOUND_HPP

#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
template <typename DistributedBoundaryType>
void Problem::distributed_boundary_send_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.ze, boundary.ze_at_gp);
    dbound.ComputeUgp(state.qx, boundary.qx_at_gp);
    dbound.ComputeUgp(state.qy, boundary.qy_at_gp);

    dbound.boundary_condition.SetEX(boundary.ze_at_gp, boundary.qx_at_gp, boundary.qy_at_gp);

    dbound.boundary_condition.SetWetDryEX(dbound.data.wet_dry_state.wet);
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    auto& wd_state_in = dbound.data.wet_dry_state;

    bool wet_ex;
    dbound.boundary_condition.GetWetDryEX(wet_ex);

    if (wd_state_in.wet || wet_ex) {
        const uint stage = stepper.GetStage();

        auto& state    = dbound.data.state[stage];
        auto& boundary = dbound.data.boundary[dbound.bound_id];
        auto& sp_at_gp = dbound.data.spherical_projection.sp_at_gp_boundary[dbound.bound_id];

        dbound.boundary_condition.ComputeFlux(stepper, dbound);

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < dbound.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] -= dbound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
            state.rhs_qx[dof] -= dbound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
            state.rhs_qy[dof] -= dbound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
        }
    }
}
}

#endif
