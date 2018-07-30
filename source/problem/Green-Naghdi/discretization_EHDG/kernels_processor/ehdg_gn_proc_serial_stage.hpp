#ifndef EHDG_GN_PROC_SERIAL_STAGE_HPP
#define EHDG_GN_PROC_SERIAL_STAGE_HPP

#include "general_definitions.hpp"

#include "ehdg_gn_kernels_processor.hpp"
#include "ehdg_gn_proc_derivatives.hpp"
namespace GN {
namespace EHDG {
void Problem::serial_stage_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    Problem::serial_swe_stage_kernel(stepper, discretization);

    Problem::serial_dispersive_correction_kernel(stepper, discretization);
}

void Problem::serial_swe_stage_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    /* Global Step */
    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::global_swe_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { Problem::global_swe_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::global_swe_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::global_swe_edge_boundary_kernel(stepper, edge_bound); });
    /* Global Step */

    /* Local Step */
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_swe_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_swe_source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::local_swe_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { Problem::local_swe_boundary_kernel(stepper, bound); });
    /* Local Step */

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::update_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = Problem::scrutinize_solution_kernel(stepper, elt);

        if (nan_found)
            abort();
    });
}

void Problem::serial_dispersive_correction_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    // Compute du, ddu
    Problem::serial_derivatives_kernel(stepper, discretization);
}
}
}

#endif