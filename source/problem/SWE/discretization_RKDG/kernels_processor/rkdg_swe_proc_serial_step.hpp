#ifndef RKDG_SWE_PROC_SERIAL_STEP_HPP
#define RKDG_SWE_PROC_SERIAL_STEP_HPP

#include "rkdg_swe_kernels_processor.hpp"
#include "problem/SWE/problem_slope_limiter/swe_CS_sl_serial.hpp"
#include "problem/SWE/seabed_update/swe_seabed_update.hpp"

namespace SWE {
namespace RKDG {
template <template <typename> class DiscretizationType, typename ProblemType>
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

template <template <typename> class DiscretizationType, typename ProblemType>
void Problem::stage_serial(DiscretizationType<ProblemType>& discretization,
                           typename ProblemType::ProblemGlobalDataType& global_data,
                           ProblemStepperType& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) { Problem::boundary_kernel(stepper, bound); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& state = elt.data.state[stepper.GetStage()];

        state.solution = elt.ApplyMinv(state.rhs);

        stepper.UpdateState(elt);
    });

    if (SWE::SedimentTransport::bed_update)
        SWE::seabed_update(stepper, discretization);

    ++stepper;

    if (SWE::SedimentTransport::bed_slope_limiting)
        SWE::CS_seabed_slope_limiter(stepper, discretization);

    if (SWE::SedimentTransport::bed_update)
        SWE::seabed_data_update(stepper, discretization);

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
