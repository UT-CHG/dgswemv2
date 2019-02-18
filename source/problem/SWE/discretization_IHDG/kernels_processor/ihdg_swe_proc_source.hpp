#ifndef IHDG_SWE_PROC_SOURCE_HPP
#define IHDG_SWE_PROC_SOURCE_HPP

#include "problem/SWE/problem_souce/swe_source.hpp"

namespace SWE {
namespace IHDG {
template <typename StepperType, typename ElementType>
void Problem::init_source_kernel(const StepperType& stepper, ElementType& elt) {
    if (stepper.GetOrder() == 2) {
        auto& internal = elt.data.internal;

        SWE::get_source(stepper.GetTimeAtCurrentStage(), elt);

        internal.rhs_prev += elt.IntegrationPhi(internal.source_at_gp);
    }
}

template <typename StepperType, typename ElementType>
void Problem::local_source_kernel(const StepperType& stepper, ElementType& elt) {
    auto& internal = elt.data.internal;

    SWE::get_source(stepper.GetTimeAtCurrentStage() + stepper.GetDT(), elt);

    internal.rhs_local += (1.0 - stepper.GetTheta()) * elt.IntegrationPhi(internal.source_at_gp);
}
}
}

#endif
