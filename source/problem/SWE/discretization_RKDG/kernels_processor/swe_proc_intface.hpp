#ifndef SWE_PROC_INTFACE_HPP
#define SWE_PROC_INTFACE_HPP

#include "problem/SWE/discretization_RKDG/numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
template <typename InterfaceType>
void Problem::interface_kernel(const RKStepper& stepper, InterfaceType& intface) {
    auto& wd_state_in = intface.data_in.wet_dry_state;
    auto& wd_state_ex = intface.data_ex.wet_dry_state;

    if (wd_state_in.wet || wd_state_ex.wet) {
        const uint stage = stepper.GetStage();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        intface.ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
        intface.ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
        intface.ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

        intface.ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
        intface.ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
        intface.ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);

        intface.specialization.ComputeFlux(stepper, intface);

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
            state_in.rhs_ze[dof] -= intface.IntegrationPhiIN(dof, boundary_in.ze_numerical_flux_at_gp);
            state_in.rhs_qx[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qx_numerical_flux_at_gp);
            state_in.rhs_qy[dof] -= intface.IntegrationPhiIN(dof, boundary_in.qy_numerical_flux_at_gp);
        }

        for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
            state_ex.rhs_ze[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.ze_numerical_flux_at_gp);
            state_ex.rhs_qx[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qx_numerical_flux_at_gp);
            state_ex.rhs_qy[dof] -= intface.IntegrationPhiEX(dof, boundary_ex.qy_numerical_flux_at_gp);
        }
    }
}
}
}

#endif
