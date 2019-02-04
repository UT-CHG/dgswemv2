#ifndef EHDG_SWE_PROC_OMPI_STEP_HPP
#define EHDG_SWE_PROC_OMPI_STEP_HPP

#include "ehdg_swe_kernels_processor.hpp"

namespace SWE {
namespace EHDG {
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
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Current (time, stage): (" << sim->stepper.GetTimeAtCurrentStage()
                                                  << ',' << sim->stepper.GetStage() << ')' << std::endl;

            sim_units[su_id]->writer.GetLogFile() << "Exchanging data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(CommTypes::bound_state, sim->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [sim](auto& dbound) { Problem::global_distributed_boundary_kernel(sim->stepper, dbound); });

        sim_units[su_id]->communicator.SendAll(CommTypes::bound_state, sim->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        /* Global Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [sim](auto& intface) { Problem::global_interface_kernel(sim->stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [sim](auto& bound) { Problem::global_boundary_kernel(sim->stepper, bound); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [sim](auto& edge_int) { Problem::global_edge_interface_kernel(sim->stepper, edge_int); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [sim](auto& edge_bound) { Problem::global_edge_boundary_kernel(sim->stepper, edge_bound); });
        /* Global Pre Receive Step */

        /* Local Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [sim](auto& elt) { Problem::local_volume_kernel(sim->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [sim](auto& elt) { Problem::local_source_kernel(sim->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [sim](auto& intface) { Problem::local_interface_kernel(sim->stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [sim](auto& bound) { Problem::local_boundary_kernel(sim->stepper, bound); });
        /* Local Pre Receive Step */

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Starting to wait on receive with timestamp: " << sim->stepper.GetTimestamp() << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::bound_state, sim->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        /* Global Post Receive Step */
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [sim](auto& edge_dbound) { Problem::global_edge_distributed_kernel(sim->stepper, edge_dbound); });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [sim](auto& dbound) { Problem::local_distributed_boundary_kernel(sim->stepper, dbound); });

        sim_units[su_id]->discretization.mesh.CallForEachElement([sim](auto& elt) {
            auto& state = elt.data.state[sim->stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            sim->stepper.UpdateState(elt);
        });
        /* Local Post Receive Step */

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

    ++(sim->stepper);

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

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::bound_state, sim->stepper.GetTimestamp());
    }
}
}
}

#endif
