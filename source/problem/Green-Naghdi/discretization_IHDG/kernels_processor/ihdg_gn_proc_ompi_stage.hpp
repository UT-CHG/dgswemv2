#ifndef IHDG_GN_PROC_OMPI_STAGE_HPP
#define IHDG_GN_PROC_OMPI_STAGE_HPP

#include "general_definitions.hpp"

#include "ihdg_gn_kernels_processor.hpp"

namespace GN {
namespace IHDG {
template <typename OMPISimUnitType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

    n_threads = (uint)omp_get_num_threads();
    thread_id = (uint)omp_get_thread_num();

    sim_per_thread = (sim_units.size() + n_threads - 1) / n_threads;

    begin_sim_id = sim_per_thread * thread_id;
    end_sim_id   = std::min(sim_per_thread * (thread_id + 1), (uint)sim_units.size());

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Current (time, stage): (" << sim_units[su_id]->stepper.GetTimeAtCurrentStage() << ','
                << sim_units[su_id]->stepper.GetStage() << ')' << std::endl;

            sim_units[su_id]->writer.GetLogFile() << "Exchanging data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::global_swe_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->communicator.SendAll(sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        if (sim_units[su_id]->parser.ParsingInput()) {
            if (sim_units[su_id]->writer.WritingVerboseLog()) {
                sim_units[su_id]->writer.GetLogFile() << "Parsing input" << std::endl;
            }

            sim_units[su_id]->parser.ParseInput(sim_units[su_id]->stepper, sim_units[su_id]->discretization.mesh);
        }

        /* Global Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            Problem::global_swe_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
            Problem::global_swe_boundary_kernel(sim_units[su_id]->stepper, bound);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&sim_units, su_id](auto& edge_int) {
            Problem::global_swe_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&sim_units, su_id](auto& edge_bound) {
            Problem::global_swe_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
        });
        /* Global Pre Receive Step */

        /* Local Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::local_swe_volume_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::local_swe_source_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            Problem::local_swe_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&sim_units, su_id](auto& bound) { Problem::local_swe_boundary_kernel(sim_units[su_id]->stepper, bound); });
        /* Local Pre Receive Step */

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Starting to wait on receive with timestamp: " << sim_units[su_id]->stepper.GetTimestamp()
                << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        /* Global Post Receive Step */
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&sim_units, su_id](auto& edge_dbound) {
                Problem::global_swe_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
            });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::local_swe_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
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

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(sim_units[su_id]->stepper.GetTimestamp());
    }
}
}
}

#endif