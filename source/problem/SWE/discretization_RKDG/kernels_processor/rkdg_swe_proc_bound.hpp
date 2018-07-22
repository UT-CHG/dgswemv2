#ifndef RKDG_SWE_PROC_BOUND_HPP
#define RKDG_SWE_PROC_BOUND_HPP

#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
template <typename BoundaryType>
void Problem::boundary_kernel(const RKStepper& stepper, BoundaryType& bound) {
    auto& wd_state = bound.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.q_at_gp = bound.ComputeUgp(state.q);

        bound.boundary_condition.ComputeFlux(
            stepper, bound.surface_normal, boundary.q_at_gp, boundary.aux_at_gp, boundary.F_hat_at_gp);

        // now compute contributions to the righthand side
        for (uint dof = 0; dof < bound.data.get_ndof(); ++dof) {
            column(state.rhs, dof) -= bound.IntegrationPhi(dof, boundary.F_hat_at_gp);
        }
    }
}
}
}

#endif
