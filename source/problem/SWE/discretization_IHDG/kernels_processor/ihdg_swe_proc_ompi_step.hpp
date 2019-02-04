#ifndef IHDG_SWE_PROC_OMPI_STEP_HPP
#define IHDG_SWE_PROC_OMPI_STEP_HPP

#include "ihdg_swe_kernels_processor.hpp"
#include "ihdg_swe_proc_init_iter.hpp"
#include "ihdg_swe_proc_ompi_sol_glob_prob.hpp"

namespace SWE {
namespace IHDG {
template <typename OMPISimType>
void Problem::step_ompi(OMPISimType* sim, uint begin_sim_id, uint end_sim_id) {
    auto& sim_units = sim->sim_units;
    // Here one assumes that there is at lease one sim unit present
    // This is of course not always true
    for (uint stage = 0; stage < sim->stepper.GetNumStages(); ++stage) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (sim_units[su_id]->parser.ParsingInput()) {
                sim_units[su_id]->parser.ParseInput(sim->stepper, sim_units[su_id]->discretization.mesh);
            }
        }

        Problem::stage_ompi(sim, begin_sim_id, end_sim_id);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingOutput()) {
            sim_units[su_id]->writer.WriteOutput(sim->stepper, sim_units[su_id]->discretization.mesh);
        }
    }
}

template <typename OMPISimType>
void Problem::stage_ompi(OMPISimType* sim, uint begin_sim_id, uint end_sim_id) {
    auto& sim_units = sim->sim_units;

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        Problem::init_iteration(sim->stepper, sim_units[su_id]->discretization);
    }

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            /* Local Step */
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [sim](auto& elt) { Problem::local_volume_kernel(sim->stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [sim](auto& elt) { Problem::local_source_kernel(sim->stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachInterface(
                [sim](auto& intface) { Problem::local_interface_kernel(sim->stepper, intface); });

            sim_units[su_id]->discretization.mesh.CallForEachBoundary(
                [sim](auto& bound) { Problem::local_boundary_kernel(sim->stepper, bound); });

            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
                [sim](auto& dbound) { Problem::local_distributed_boundary_kernel(sim->stepper, dbound); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [sim](auto& edge_int) { Problem::local_edge_interface_kernel(sim->stepper, edge_int); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [sim](auto& edge_bound) { Problem::local_edge_boundary_kernel(sim->stepper, edge_bound); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [sim](auto& edge_dbound) { Problem::local_edge_distributed_kernel(sim->stepper, edge_dbound); });
            /* Local Step */

            /* Global Step */
            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [sim](auto& edge_int) { Problem::global_edge_interface_kernel(sim->stepper, edge_int); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [sim](auto& edge_bound) { Problem::global_edge_boundary_kernel(sim->stepper, edge_bound); });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [sim](auto& edge_dbound) { Problem::global_edge_distributed_kernel(sim->stepper, edge_dbound); });
            /* Global Step */
        }

        bool converged = ompi_solve_global_problem(sim, begin_sim_id, end_sim_id);

        if (converged) {
            break;
        }
    }

    ++(sim->stepper);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([sim](auto& elt) {
            uint n_stages = sim->stepper.GetNumStages();

            auto& state = elt.data.state;

            std::swap(state[0].q, state[n_stages].q);
        });
    }

    if (SWE::PostProcessing::slope_limiting) {
        CS_slope_limiter_ompi(sim->stepper, sim->sim_units, begin_sim_id, end_sim_id, CommTypes::baryctr_state);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([sim](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(sim->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });
    }
}
}
}

#endif
