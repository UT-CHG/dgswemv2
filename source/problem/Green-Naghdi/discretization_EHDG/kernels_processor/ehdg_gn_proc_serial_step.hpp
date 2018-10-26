#ifndef EHDG_GN_PROC_SERIAL_STEP_HPP
#define EHDG_GN_PROC_SERIAL_STEP_HPP

#include "general_definitions.hpp"

#include "ehdg_gn_kernels_processor.hpp"
#include "ehdg_gn_proc_serial_sol_glob_prob.hpp"
#include "ehdg_gn_proc_derivatives_serial.hpp"

namespace GN {
namespace EHDG {
template <typename SerialSimType>
void Problem::step_serial(SerialSimType* sim) {
    for (uint stage = 0; stage < sim->stepper.GetNumStages(); ++stage) {
        if (sim->parser.ParsingInput()) {
            sim->parser.ParseInput(sim->stepper, sim->discretization.mesh);
        }

        Problem::stage_serial(sim->stepper, sim->discretization);

        ++(sim->stepper);
    }

    sim->discretization.mesh.CallForEachElement([sim](auto& elt) { Problem::swap_states_kernel(sim->stepper, elt); });

    if (sim->writer.WritingOutput()) {
        sim->writer.WriteOutput(sim->stepper, sim->discretization.mesh);
    }
}

void Problem::stage_serial(const ProblemStepperType& stepper, ProblemDiscretizationType& discretization) {
    Problem::swe_stage_serial(stepper, discretization);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& curr_state = elt.data.state[stage];
        auto& next_state = elt.data.state[stage + 1];

        curr_state.q = next_state.q;
    });

    Problem::dispersive_correction_serial(stepper, discretization);

    Problem::swe_stage_serial(stepper, discretization);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = SWE::EHDG::Problem::scrutinize_solution_kernel(stepper, elt);

        if (nan_found)
            abort();
    });
}

void Problem::swe_stage_serial(const ProblemStepperType& stepper, ProblemDiscretizationType& discretization) {
    /* Global Step */
    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { SWE::EHDG::Problem::global_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { SWE::EHDG::Problem::global_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { SWE::EHDG::Problem::global_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { SWE::EHDG::Problem::global_edge_boundary_kernel(stepper, edge_bound); });
    /* Global Step */

    /* Local Step */
    discretization.mesh.CallForEachElement(
        [&stepper](auto& elt) { SWE::EHDG::Problem::local_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement(
        [&stepper](auto& elt) { SWE::EHDG::Problem::local_source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { SWE::EHDG::Problem::local_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary(
        [&stepper](auto& bound) { SWE::EHDG::Problem::local_boundary_kernel(stepper, bound); });
    /* Local Step */

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& curr_state = elt.data.state[stage];
        auto& next_state = elt.data.state[stage + 1];

        curr_state.solution = elt.ApplyMinv(curr_state.rhs);

        next_state.q = curr_state.q + stepper.GetDT() / 2.0 * curr_state.solution;
    });
}

void Problem::dispersive_correction_serial(const ProblemStepperType& stepper,
                                           ProblemDiscretizationType& discretization) {
    // Compute du, ddu
    Problem::compute_derivatives_serial(stepper, discretization);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_dc_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_dc_source_kernel(stepper, elt); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::local_dc_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::local_dc_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::global_dc_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::global_dc_edge_boundary_kernel(stepper, edge_bound); });

    Problem::serial_solve_global_dc_problem(stepper, discretization);

    discretization.mesh.CallForEachElement(
        [&stepper](auto& elt) { Problem::dispersive_correction_kernel(stepper, elt); });
}
}
}

#endif