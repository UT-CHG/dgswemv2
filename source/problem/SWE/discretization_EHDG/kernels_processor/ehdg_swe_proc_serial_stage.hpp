#ifndef EHDG_SWE_PROC_SERIAL_STAGE_HPP
#define EHDG_SWE_PROC_SERIAL_STAGE_HPP

#include "general_definitions.hpp"

#include "ehdg_swe_kernels_processor.hpp"

namespace SWE {
namespace EHDG {
void Problem::stage_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    /* Global Step */
    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::global_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { Problem::global_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::global_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    /* Global Step */

    /* Local Step */
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::local_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { Problem::local_boundary_kernel(stepper, bound); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::update_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = Problem::scrutinize_solution_kernel(stepper, elt);

        if (nan_found)
            abort();
    });
    /* Local Step */
}
}
}

#endif