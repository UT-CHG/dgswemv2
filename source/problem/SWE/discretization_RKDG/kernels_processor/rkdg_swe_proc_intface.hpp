#ifndef RKDG_SWE_PROC_INTFACE_HPP
#define RKDG_SWE_PROC_INTFACE_HPP

#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

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

        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);

        intface.specialization.ComputeFlux(stepper, intface);

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < intface.data_in.get_ndof(); ++dof) {
            column(state_in.rhs, dof) -= intface.IntegrationPhiIN(dof, boundary_in.F_hat_at_gp);
        }

        for (uint dof = 0; dof < intface.data_ex.get_ndof(); ++dof) {
            column(state_ex.rhs, dof) -= intface.IntegrationPhiEX(dof, boundary_ex.F_hat_at_gp);
        }
    }
}
}
}

#endif
