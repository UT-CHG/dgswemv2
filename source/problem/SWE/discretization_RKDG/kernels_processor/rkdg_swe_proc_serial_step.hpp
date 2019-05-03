#ifndef RKDG_SWE_PROC_SERIAL_STEP_HPP
#define RKDG_SWE_PROC_SERIAL_STEP_HPP

#include "rkdg_swe_kernels_processor.hpp"
#include "problem/SWE/problem_slope_limiter/swe_CS_sl_serial.hpp"

namespace SWE {
namespace RKDG {
template <template <typename> typename DiscretizationType, typename ProblemType>
void Problem::step_serial(DiscretizationType<ProblemType>& discretization,
                          typename ProblemType::ProblemGlobalDataType& global_data,
                          ProblemStepperType& stepper,
                          typename ProblemType::ProblemWriterType& writer,
                          typename ProblemType::ProblemParserType& parser) {
    for (uint stage = 0; stage < stepper.GetNumStages(); ++stage) {
        if (parser.ParsingInput()) {
            parser.ParseInput(stepper, discretization.mesh);
        }

        Problem::stage_serial(discretization, global_data, stepper);
    }

    if (writer.WritingOutput()) {
        writer.WriteOutput(stepper, discretization.mesh);
    }
}

template <template <typename> typename DiscretizationType, typename ProblemType>
void Problem::stage_serial(DiscretizationType<ProblemType>& discretization,
                           typename ProblemType::ProblemGlobalDataType& global_data,
                           ProblemStepperType& stepper) {
    VolumeKernel<ProblemStepperType> volume_kernel(stepper);
    discretization.mesh.CallForEachElement(volume_kernel);

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::source_kernel(stepper, elt); });

    InterfaceKernel interface_kernel(stepper);
    discretization.mesh.CallForEachInterface(interface_kernel);

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) { Problem::boundary_kernel(stepper, bound); });

    UpdateKernel update_kernel(stepper);
    discretization.mesh.CallForEachElement(update_kernel);

    ++stepper;

    if (SWE::PostProcessing::wetting_drying) {
        discretization.mesh.CallForEachElement([&stepper](auto& elt) { wetting_drying_kernel(stepper, elt); });
    }

    if (SWE::PostProcessing::slope_limiting) {
        CS_slope_limiter_serial(stepper, discretization);
    }

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = SWE::scrutinize_solution(stepper, elt);

        if (nan_found) {
            std::cerr << "Fatal Error: NaN found at element " << elt.GetID() << std::endl;
            abort();
        }
    });
}
}
}

#endif