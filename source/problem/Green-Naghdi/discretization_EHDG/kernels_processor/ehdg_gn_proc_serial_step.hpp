#ifndef EHDG_GN_PROC_SERIAL_STEP_HPP
#define EHDG_GN_PROC_SERIAL_STEP_HPP

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

        SWE_SIM::Problem::stage_serial(sim->stepper, sim->discretization);
    }

    for (uint stage = 0; stage < sim->stepper.GetNumStages(); ++stage) {
        if (sim->parser.ParsingInput()) {
            sim->parser.ParseInput(sim->stepper, sim->discretization.mesh);
        }

        Problem::dispersive_correction_serial(sim->stepper, sim->discretization);
    }

    for (uint stage = 0; stage < sim->stepper.GetNumStages(); ++stage) {
        if (sim->parser.ParsingInput()) {
            sim->parser.ParseInput(sim->stepper, sim->discretization.mesh);
        }

        SWE_SIM::Problem::stage_serial(sim->stepper, sim->discretization);
    }

    if (sim->writer.WritingOutput()) {
        sim->writer.WriteOutput(sim->stepper, sim->discretization.mesh);
    }
}

void Problem::dispersive_correction_serial(ProblemStepperType& stepper,
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

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        Problem::dispersive_correction_kernel(stepper, elt);

        auto& state = elt.data.state[stepper.GetStage()];

        state.solution = elt.ApplyMinv(state.rhs);

        stepper.UpdateState(elt);
    });

    ++stepper;
}
}
}

#endif