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

    boundary.q_at_gp = dbound.ComputeUgp(state.q);

    // Construct message to exterior state
    std::vector<double> message;

    message.reserve(1 + SWE::n_variables * dbound.data.get_ngp_boundary(dbound.bound_id));

    message.push_back(dbound.data.wet_dry_state.wet);

    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        for (uint var = 0; var < SWE::n_variables; ++var) {
            message.push_back(boundary.q_at_gp(var, gp));
        }
    }

    // Set message to send buffer
    dbound.boundary_condition.exchanger.SetToSendBuffer(SWE::CommTypes::processor, message);
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
        state.rhs -= dbound.IntegrationPhi(boundary.F_hat_at_gp);
    }
}
}
}

#endif
