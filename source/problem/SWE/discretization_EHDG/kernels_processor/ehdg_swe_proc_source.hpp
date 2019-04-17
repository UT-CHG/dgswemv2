#ifndef EHDG_SWE_PROC_SOURCE_HPP
#define EHDG_SWE_PROC_SOURCE_HPP

#include "problem/SWE/problem_source/swe_source.hpp"

namespace SWE {
namespace EHDG {
template <typename ElementType>
void Problem::local_source_kernel(const ProblemStepperType& stepper, ElementType& elt) {
    if (elt.data.wet_dry_state.wet) {
        auto& state    = elt.data.state[stepper.GetStage()];
        auto& internal = elt.data.internal;

        SWE::get_source(stepper.GetTimeAtCurrentStage(), elt);

        state.rhs += elt.IntegrationPhi(internal.source_at_gp);
    }
}
}
}

#endif
