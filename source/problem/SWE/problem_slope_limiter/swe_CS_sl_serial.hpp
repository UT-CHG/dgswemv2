#ifndef SWE_CS_SL_SERIAL_HPP
#define SWE_CS_SL_SERIAL_HPP

// Implementation of Cockburn-Shu slope limiter
#include "swe_CS_slope_limiter.hpp"

namespace SWE {
template <typename StepperType, typename DiscretizationType>
void CS_slope_limiter_serial(StepperType& stepper, DiscretizationType& discretization) {
    discretization.mesh.CallForEachElement(
        [&stepper](auto& elt) { slope_limiting_prepare_element_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { slope_limiting_prepare_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { slope_limiting_prepare_boundary_kernel(stepper, bound); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { slope_limiting_kernel(stepper, elt); });
}
}

#endif