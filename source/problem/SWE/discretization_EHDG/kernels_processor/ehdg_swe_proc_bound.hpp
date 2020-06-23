#ifndef EHDG_SWE_PROC_BOUND_HPP
#define EHDG_SWE_PROC_BOUND_HPP

namespace SWE {
namespace EHDG {
template <typename BoundaryType>
void Problem::global_boundary_kernel(const ProblemStepperType& stepper, BoundaryType& bound) {
    if (bound.data.wet_dry_state.wet) {
        auto& state      = bound.data.state[stepper.GetStage()];
        auto& boundary   = bound.data.boundary[bound.bound_id];
        boundary.q_at_gp = bound.ComputeUgp(state.q);
        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);
    }
}

template <typename BoundaryType>
void Problem::local_boundary_kernel(const ProblemStepperType& stepper, BoundaryType& bound) {
    if (bound.data.wet_dry_state.wet) {
        auto& state    = bound.data.state[stepper.GetStage()];
        auto& boundary = bound.data.boundary[bound.bound_id];
        state.rhs -= bound.IntegrationPhi(boundary.F_hat_at_gp);
    }
}
}
}

#endif
