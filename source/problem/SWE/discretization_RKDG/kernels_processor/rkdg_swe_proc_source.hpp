#ifndef RKDG_SWE_PROC_SOURCE_HPP
#define RKDG_SWE_PROC_SOURCE_HPP

#include "problem/SWE/problem_source/swe_source.hpp"

namespace SWE {
namespace RKDG {
template <typename StepperType, typename ElementType>
void Problem::source_kernel(const StepperType& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        SWE::get_source(stepper.GetTimeAtCurrentStage(), elt);

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            state.rhs[var] += elt.IntegrationPhi(row(internal.source_at_gp,var));
        }
    }
}
}
}

#endif
