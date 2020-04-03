#ifndef EHDG_GN_PROC_SERIAL_STEP_HPP
#define EHDG_GN_PROC_SERIAL_STEP_HPP

#include "ehdg_gn_kernels_processor.hpp"
#include "ehdg_gn_proc_serial_sol_glob_prob.hpp"
#include "ehdg_gn_proc_derivatives_serial.hpp"

namespace GN {
namespace EHDG {
template <template <typename> class DiscretizationType, typename ProblemType>
void seabed_update(DiscretizationType<ProblemType>& discretization,
                   typename ProblemType::ProblemGlobalDataType& global_data,
                   typename ProblemType::ProblemStepperType& stepper);

void Problem::step_serial(ProblemDiscretizationType& discretization,
                          ProblemGlobalDataType& global_data,
                          ProblemStepperType& stepper,
                          ProblemWriterType& writer,
                          ProblemParserType& parser) {
    auto& first_stepper  = stepper.GetFirstStepper();
    auto& second_stepper = stepper.GetSecondStepper();

    for (uint stage = 0; stage < first_stepper.GetNumStages(); ++stage) {
        if (parser.ParsingInput()) {
            parser.ParseInput(first_stepper, discretization.mesh);
        }

        SWE_SIM::Problem::stage_serial(discretization, global_data, first_stepper);
    }

    for (uint stage = 0; stage < second_stepper.GetNumStages(); ++stage) {
        if (parser.ParsingInput()) {
            parser.ParseInput(second_stepper, discretization.mesh);
        }

        Problem::dispersive_correction_serial(discretization, global_data, second_stepper);
    }

    for (uint stage = 0; stage < first_stepper.GetNumStages(); ++stage) {
        if (parser.ParsingInput()) {
            parser.ParseInput(first_stepper, discretization.mesh);
        }

        SWE_SIM::Problem::stage_serial(discretization, global_data, first_stepper);
    }

    if (SWE::PostProcessing::bed_update && (second_stepper.GetStep() % SWE::PostProcessing::bed_update_freq) == 0) {
        seabed_update(discretization, global_data, stepper);
    }

    if (writer.WritingOutput()) {
        writer.WriteOutput(stepper, discretization.mesh);
    }
}

void Problem::dispersive_correction_serial(ProblemDiscretizationType& discretization,
                                           ProblemGlobalDataType& global_data,
                                           ESSPRKStepper& stepper) {
    // Compute du, ddu
    Problem::compute_derivatives_serial(discretization, global_data, stepper);

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

    Problem::serial_solve_global_dc_problem(discretization, global_data, stepper);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        Problem::dispersive_correction_kernel(stepper, elt);

        auto& state = elt.data.state[stepper.GetStage()];

        state.solution = elt.ApplyMinv(state.rhs);

        stepper.UpdateState(elt);
    });

    ++stepper;
}

template <template <typename> class DiscretizationType, typename ProblemType>
void seabed_update(DiscretizationType<ProblemType>& discretization,
                   typename ProblemType::ProblemGlobalDataType& global_data,
                   typename ProblemType::ProblemStepperType& stepper) {
    SWE::seabed_update(stepper.GetSecondStepper(), discretization);
    if (SWE::PostProcessing::bed_slope_limiting)
        SWE::CS_seabed_slope_limiter(discretization);
    SWE::seabed_data_update(discretization);
    ProblemType::compute_bathymetry_derivatives_serial(discretization, global_data);
}
}
}

#endif