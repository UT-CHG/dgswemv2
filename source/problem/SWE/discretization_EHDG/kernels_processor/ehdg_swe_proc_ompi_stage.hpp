#ifndef EHDG_SWE_PROC_OMPI_STAGE_HPP
#define EHDG_SWE_PROC_OMPI_STAGE_HPP

#include "general_definitions.hpp"

#include "ehdg_swe_kernels_processor.hpp"

namespace SWE {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, uint begin_sim_id, uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Current (time, stage): (" << sim_units[su_id]->stepper.GetTimeAtCurrentStage() << ','
                << sim_units[su_id]->stepper.GetStage() << ')' << std::endl;

            sim_units[su_id]->writer.GetLogFile() << "Exchanging data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(CommTypes::bound_state, sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::global_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::bound_state, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        /* Global Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            Problem::global_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&sim_units, su_id](auto& bound) { Problem::global_boundary_kernel(sim_units[su_id]->stepper, bound); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&sim_units, su_id](auto& edge_int) {
            Problem::global_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&sim_units, su_id](auto& edge_bound) {
            Problem::global_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
        });
        /* Global Pre Receive Step */

        /* Local Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::local_volume_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::local_source_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            Problem::local_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&sim_units, su_id](auto& bound) { Problem::local_boundary_kernel(sim_units[su_id]->stepper, bound); });
        /* Local Pre Receive Step */

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Starting to wait on receive with timestamp: " << sim_units[su_id]->stepper.GetTimestamp()
                << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::bound_state,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        /* Global Post Receive Step */
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&sim_units, su_id](auto& edge_dbound) {
                Problem::global_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
            });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::local_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::update_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = Problem::scrutinize_solution_kernel(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });
        /* Local Post Receive Step */

        ++(sim_units[su_id]->stepper);

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::bound_state, sim_units[su_id]->stepper.GetTimestamp());
    }
}
}
}

#endif