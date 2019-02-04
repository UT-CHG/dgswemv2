#ifndef IHDG_SWE_PROC_OMPI_STEP_HPP
#define IHDG_SWE_PROC_OMPI_STEP_HPP

#include "ihdg_swe_kernels_processor.hpp"
#include "ihdg_swe_proc_init_iter.hpp"
#include "ihdg_swe_proc_ompi_sol_glob_prob.hpp"

namespace SWE {
namespace IHDG {
template <template <typename> typename OMPISimUnitType, typename ProblemType>
void Problem::step_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                        typename ProblemType::ProblemGlobalDataType& global_data,
                        typename ProblemType::ProblemStepperType& stepper,
                        uint begin_sim_id,
                        uint end_sim_id) {
    for (uint stage = 0; stage < stepper.GetNumStages(); ++stage) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (sim_units[su_id]->parser.ParsingInput()) {
                sim_units[su_id]->parser.ParseInput(stepper, sim_units[su_id]->discretization.mesh);
            }
        }

        Problem::stage_ompi(sim_units, global_data, stepper, begin_sim_id, end_sim_id);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingOutput()) {
            sim_units[su_id]->writer.WriteOutput(stepper, sim_units[su_id]->discretization.mesh);
        }
    }
}

template <template <typename> typename OMPISimUnitType, typename ProblemType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                         typename ProblemType::ProblemGlobalDataType& global_data,
                         typename ProblemType::ProblemStepperType& stepper,
                         uint begin_sim_id,
                         uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        Problem::init_iteration(stepper, sim_units[su_id]->discretization);
    }

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            /* Local Step */
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&stepper](auto& elt) { Problem::local_volume_kernel(stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&stepper](auto& elt) { Problem::local_source_kernel(stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachInterface(
                [&stepper](auto& intface) { Problem::local_interface_kernel(stepper, intface); });

            sim_units[su_id]->discretization.mesh.CallForEachBoundary(
                [&stepper](auto& bound) { Problem::local_boundary_kernel(stepper, bound); });

            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
                [&stepper](auto& dbound) { Problem::local_distributed_boundary_kernel(stepper, dbound); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&stepper](auto& edge_int) { Problem::local_edge_interface_kernel(stepper, edge_int); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&stepper](auto& edge_bound) { Problem::local_edge_boundary_kernel(stepper, edge_bound); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&stepper](auto& edge_dbound) { Problem::local_edge_distributed_kernel(stepper, edge_dbound); });
            /* Local Step */

            /* Global Step */
            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&stepper](auto& edge_int) { Problem::global_edge_interface_kernel(stepper, edge_int); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&stepper](auto& edge_bound) { Problem::global_edge_boundary_kernel(stepper, edge_bound); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&stepper](auto& edge_dbound) { Problem::global_edge_distributed_kernel(stepper, edge_dbound); });
            /* Global Step */
        }

        bool converged = ompi_solve_global_problem(sim_units, global_data, stepper, begin_sim_id, end_sim_id);

        if (converged) {
            break;
        }
    }

#pragma omp barrier
#pragma omp master
    {
    ++(stepper);
    }
#pragma omp barrier

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            uint n_stages = stepper.GetNumStages();

            auto& state = elt.data.state;

            std::swap(state[0].q, state[n_stages].q);
        });
    }

    if (SWE::PostProcessing::slope_limiting) {
        CS_slope_limiter_ompi(stepper, sim_units, begin_sim_id, end_sim_id, CommTypes::baryctr_state);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });
    }
}
}
}

#endif
