#ifndef RKDG_SWE_PROC_VOLUME_HPP
#define RKDG_SWE_PROC_VOLUME_HPP

#include "problem/SWE/problem_function_files/swe_compute_F.hpp"

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::volume_kernel(const RKStepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);

        SWE::compute_F(internal.q_at_gp, internal.aux_at_gp, internal.Fx_at_gp, internal.Fy_at_gp);

        state.rhs = elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) +
                    elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
    }
}
}
}

#endif
