#ifndef SWE_PROC_DBOUND_HPP
#define SWE_PROC_DBOUND_HPP

#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
template <typename DistributedBoundaryType>
void Problem::distributed_boundary_send_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.get_stage();

    auto& state = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.ze, boundary.ze_at_gp);
    dbound.ComputeUgp(state.qx, boundary.qx_at_gp);
    dbound.ComputeUgp(state.qy, boundary.qy_at_gp);

    dbound.boundary_condition.SetEX(boundary.ze_at_gp, boundary.qx_at_gp, boundary.qy_at_gp);
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_kernel(const Stepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.get_stage();

    auto& state = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        dbound.boundary_condition.GetEX(stepper,
                                        gp,
                                        dbound.surface_normal,
                                        boundary.ze_at_gp,
                                        boundary.qx_at_gp,
                                        boundary.qy_at_gp,
                                        ze_ex,
                                        qx_ex,
                                        qy_ex);

        LLF_flux(boundary.ze_at_gp[gp],
                 ze_ex,
                 boundary.qx_at_gp[gp],
                 qx_ex,
                 boundary.qy_at_gp[gp],
                 qy_ex,
                 boundary.bath_at_gp[gp],
                 dbound.surface_normal[gp],
                 boundary.ze_numerical_flux_at_gp[gp],
                 boundary.qx_numerical_flux_at_gp[gp],
                 boundary.qy_numerical_flux_at_gp[gp]);
    }

    // now compute contributions to the righthand side
    for (uint dof = 0; dof < dbound.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] -= dbound.IntegrationPhi(dof, boundary.ze_numerical_flux_at_gp);
        state.rhs_qx[dof] -= dbound.IntegrationPhi(dof, boundary.qx_numerical_flux_at_gp);
        state.rhs_qy[dof] -= dbound.IntegrationPhi(dof, boundary.qy_numerical_flux_at_gp);
    }
}
}

#endif
