#ifndef RKDG_SWE_PROC_BOUND_HPP
#define RKDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace RKDG {
template <typename BoundaryType>
void Problem::boundary_kernel(const ProblemStepperType& stepper, BoundaryType& bound) {
    if (bound.data.wet_dry_state.wet) {
        auto& state    = bound.data.state[stepper.GetStage()];
        auto& boundary = bound.data.boundary[bound.bound_id];

        for ( uint var = 0; var < SWE::n_variables; ++var) {
            boundary.q_at_gp[var] = bound.ComputeUgp(state.q[var]);
        }

        bound.boundary_condition.ComputeFlux(stepper, bound);

        // now compute contributions to the righthand side
        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            state.rhs[var] -= bound.IntegrationPhi(boundary.F_hat_at_gp[var]);
        }
    }
}
}
}

#endif
