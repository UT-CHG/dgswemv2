#ifndef EHDG_SWE_PROC_BOUND_HPP
#define EHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace EHDG {
template <typename StepperType, typename BoundaryType>
void Problem::global_boundary_kernel(const StepperType& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    boundary.q_at_gp = bound.ComputeUgp(state.q);
}

template <typename StepperType, typename BoundaryType>
void Problem::local_boundary_kernel(const StepperType& stepper, BoundaryType& bound) {
    const uint stage = stepper.GetStage();

    auto& state    = bound.data.state[stage];
    auto& boundary = bound.data.boundary[bound.bound_id];

    // now compute contributions to the righthand side
    state.rhs -= bound.IntegrationPhi(boundary.F_hat_at_gp);
}
}
}

#endif
