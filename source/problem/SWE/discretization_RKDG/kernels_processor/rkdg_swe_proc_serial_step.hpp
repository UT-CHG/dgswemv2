#ifndef RKDG_SWE_PROC_SERIAL_STEP_HPP
#define RKDG_SWE_PROC_SERIAL_STEP_HPP

#include "general_definitions.hpp"
#include "rkdg_swe_kernels_processor.hpp"

namespace SWE {
namespace RKDG {
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

void Problem::stage_serial(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) { Problem::boundary_kernel(stepper, bound); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::update_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = Problem::scrutinize_solution_kernel(stepper, elt);

        if (nan_found) {
            std::cerr << "Fatal Error: NaN found at element " << elt.GetID() << std::endl;
            abort();
        }
    });

    if (SWE::PostProcessing::wetting_drying) {
        discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::wetting_drying_kernel(stepper, elt); });
    }

    if (SWE::PostProcessing::slope_limiting) {
        discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::slope_limiting_prepare_element_kernel(stepper, elt); });

        discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { Problem::slope_limiting_prepare_interface_kernel(stepper, intface); });

        discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { Problem::slope_limiting_prepare_boundary_kernel(stepper, bound); });

        discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::slope_limiting_kernel(stepper, elt); });
    }
}
}
}

#endif