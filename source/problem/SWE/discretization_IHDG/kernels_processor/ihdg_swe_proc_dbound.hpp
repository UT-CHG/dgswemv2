#ifndef IHDG_SWE_PROC_DBOUND_HPP
#define IHDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename DistributedBoundaryType>
void Problem::init_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound) {
    if (stepper.GetOrder() == 2) {
        const uint stage = stepper.GetStage();

        auto& state_prev = dbound.data.state[stage];
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        boundary.q_at_gp = dbound.ComputeUgp(state_prev.q);
    }
}

template <typename StepperType, typename DistributedBoundaryType>
void Problem::local_distributed_boundary_kernel(const StepperType& stepper, DistributedBoundaryType& dbound) {
    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage + 1];
    auto& boundary = dbound.data.boundary[dbound.bound_id];

    boundary.q_at_gp = dbound.ComputeUgp(state.q);
}
}
}

#endif
