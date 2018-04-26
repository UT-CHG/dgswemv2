#ifndef SWE_POST_POSTPROC_HPP
#define SWE_POST_POSTPROC_HPP

#include "../../../general_definitions.hpp"

#include "swe_post_wet_dry.hpp"
#include "swe_post_slope_limit.hpp"

namespace SWE {
void Problem::postprocessor_serial_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    if (SWE::PostProcessing::slope_limiting) {
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

        mesh.CallForEachElement(slope_limiting_prepare_element_kernel);

        mesh.CallForEachInterface(slope_limiting_prepare_interface_kernel);

        mesh.CallForEachBoundary(slope_limiting_prepare_boundary_kernel);

        mesh.CallForEachElement(slope_limiting_kernel);
    }

    if (SWE::PostProcessing::wetting_drying) {
        auto wetting_drying_kernel = [&stepper](auto& elt) { Problem::wetting_drying_kernel(stepper, elt); };

        mesh.CallForEachElement(wetting_drying_kernel);
    }
}

void Problem::postprocessor_parallel_pre_send_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    if (SWE::PostProcessing::slope_limiting) {
        auto slope_limiting_prepare_element_kernel = [&stepper](auto& elt) {
            Problem::slope_limiting_prepare_element_kernel(stepper, elt);
        };

        auto slope_limiting_distributed_boundary_send_kernel = [&stepper](auto& dbound) {
            Problem::slope_limiting_distributed_boundary_send_kernel(stepper, dbound);
        };

        mesh.CallForEachElement(slope_limiting_prepare_element_kernel);

        mesh.CallForEachDistributedBoundary(slope_limiting_distributed_boundary_send_kernel);
    }
}

void Problem::postprocessor_parallel_pre_receive_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    if (SWE::PostProcessing::slope_limiting) {
        auto slope_limiting_prepare_interface_kernel = [&stepper](auto& intface) {
            Problem::slope_limiting_prepare_interface_kernel(stepper, intface);
        };

        auto slope_limiting_prepare_boundary_kernel = [&stepper](auto& bound) {
            Problem::slope_limiting_prepare_boundary_kernel(stepper, bound);
        };

        mesh.CallForEachInterface(slope_limiting_prepare_interface_kernel);

        mesh.CallForEachBoundary(slope_limiting_prepare_boundary_kernel);
    }
}

void Problem::postprocessor_parallel_post_receive_kernel(const Stepper& stepper, ProblemMeshType& mesh) {
    if (SWE::PostProcessing::slope_limiting) {
        auto slope_limiting_prepare_distributed_boundary_kernel = [&stepper](auto& dbound) {
            Problem::slope_limiting_prepare_distributed_boundary_kernel(stepper, dbound);
        };

        auto slope_limiting_kernel = [&stepper](auto& elt) { Problem::slope_limiting_kernel(stepper, elt); };

        mesh.CallForEachDistributedBoundary(slope_limiting_prepare_distributed_boundary_kernel);

        mesh.CallForEachElement(slope_limiting_kernel);
    }

    if (SWE::PostProcessing::wetting_drying) {
        auto wetting_drying_kernel = [&stepper](auto& elt) { Problem::wetting_drying_kernel(stepper, elt); };

        mesh.CallForEachElement(wetting_drying_kernel);
    }
}
}

#endif