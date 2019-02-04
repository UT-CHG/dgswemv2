#ifndef IHDG_SWE_PROC_SERIAL_STEP_HPP
#define IHDG_SWE_PROC_SERIAL_STEP_HPP

#include "ihdg_swe_kernels_processor.hpp"
#include "ihdg_swe_proc_init_iter.hpp"
#include "ihdg_swe_proc_serial_sol_glob_prob.hpp"

namespace SWE {
namespace IHDG {
template <typename ProblemType>
void Problem::step_serial(HDGDiscretization<ProblemType>& discretization,
                          typename ProblemType::ProblemGlobalDataType& global_data,
                          typename ProblemType::ProblemStepperType& stepper,
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

template <typename ProblemType>
void Problem::stage_serial(HDGDiscretization<ProblemType>& discretization,
                           typename ProblemType::ProblemGlobalDataType& global_data,
                           typename ProblemType::ProblemStepperType& stepper) {
    Problem::init_iteration(stepper, discretization);

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        /* Local Step */
        discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_volume_kernel(stepper, elt); });

        discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::local_source_kernel(stepper, elt); });

        discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { Problem::local_interface_kernel(stepper, intface); });

        discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { Problem::local_boundary_kernel(stepper, bound); });

        discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&stepper](auto& edge_int) { Problem::local_edge_interface_kernel(stepper, edge_int); });

        discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&stepper](auto& edge_bound) { Problem::local_edge_boundary_kernel(stepper, edge_bound); });

        /* Local Step */

        /* Global Step */
        discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&stepper](auto& edge_int) { Problem::global_edge_interface_kernel(stepper, edge_int); });

        discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&stepper](auto& edge_bound) { Problem::global_edge_boundary_kernel(stepper, edge_bound); });
        /* Global Step */

        bool converged = Problem::serial_solve_global_problem(discretization, global_data, stepper);

        if (converged) {
            break;
        }
    }

    ++stepper;

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        uint n_stages = stepper.GetNumStages();

        auto& state = elt.data.state;

        std::swap(state[0].q, state[n_stages].q);
    });

    if (SWE::PostProcessing::slope_limiting) {
        CS_slope_limiter_serial(stepper, discretization);
    }

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        bool nan_found = SWE::scrutinize_solution(stepper, elt);

        if (nan_found)
            abort();
    });
}
}
}

#endif