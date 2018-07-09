#ifndef RKDG_SWE_PROC_DBOUND_HPP
#define RKDG_SWE_PROC_DBOUND_HPP

#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
template <typename DistributedBoundaryType>
void Problem::distributed_boundary_send_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    dbound.ComputeUgp(state.q, boundary.q_at_gp);

    dbound.boundary_condition.exchanger.SetEX(boundary.q_at_gp);

    dbound.boundary_condition.exchanger.SetWetDryEX(dbound.data.wet_dry_state.wet);
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_kernel(const RKStepper& stepper, DistributedBoundaryType& dbound) {
    auto& wd_state_in = dbound.data.wet_dry_state;

    bool wet_ex;
    dbound.boundary_condition.exchanger.GetWetDryEX(wet_ex);

    if (wd_state_in.wet || wet_ex) {
        const uint stage = stepper.GetStage();

        auto& state    = dbound.data.state[stage];
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        dbound.boundary_condition.ComputeFlux(stepper, dbound);

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < dbound.data.get_ndof(); ++dof) {
            state.rhs[dof] -= dbound.IntegrationPhi(dof, boundary.F_hat_at_gp);
        }
    }
}
}
}

#endif
