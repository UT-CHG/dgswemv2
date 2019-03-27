#ifndef RKDG_SWE_PROC_BOUND_HPP
#define RKDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace RKDG {
template <typename StepperType, typename BoundaryType>
void Problem::boundary_kernel(const StepperType& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.q_at_gp = bound.ComputeUgp(state.q);

        bound.boundary_condition.ComputeFlux(stepper, bound);

        // now compute contributions to the righthand side
        state.rhs -= bound.IntegrationPhi(boundary.F_hat_at_gp);
    }
}
}
}

#endif
