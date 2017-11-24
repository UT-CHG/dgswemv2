#ifndef SWE_POST_STEP_POSTPROC_HPP
#define SWE_POST_STEP_POSTPROC_HPP

#include "../../../general_definitions.hpp"

#include "swe_post_wet_dry.hpp"
#include "swe_post_slope_limit.hpp"

namespace SWE {
void Problem::step_postprocessor_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    auto wetting_drying_kernel = [&stepper](auto& elt) { Problem::wetting_drying_kernel(stepper, elt); };

    auto slope_limiting_prepare_element_kernel = [&stepper](auto& elt) {
        Problem::slope_limiting_prepare_element_kernel(stepper, elt);
    };

    auto slope_limiting_prepare_interface_kernel = [&stepper](auto& intface) {
        Problem::slope_limiting_prepare_interface_kernel(stepper, intface);
    };

    auto slope_limiting_prepare_boundary_kernel = [&stepper](auto& bound) {
        Problem::slope_limiting_prepare_boundary_kernel(stepper, bound);
    };

    auto slope_limiting_kernel = [&stepper](auto& elt) { Problem::slope_limiting_kernel(stepper, elt); };

    mesh.CallForEachElement(wetting_drying_kernel);

    mesh.CallForEachElement(slope_limiting_prepare_element_kernel);

    mesh.CallForEachInterface(slope_limiting_prepare_interface_kernel);

    mesh.CallForEachBoundary(slope_limiting_prepare_boundary_kernel);

    mesh.CallForEachElement(slope_limiting_kernel);
}
}

#endif