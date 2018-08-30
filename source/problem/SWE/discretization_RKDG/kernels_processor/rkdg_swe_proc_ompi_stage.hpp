#ifndef RKDG_SWE_PROC_OMPI_STAGE_HPP
#define RKDG_SWE_PROC_OMPI_STAGE_HPP

#include "general_definitions.hpp"

#include "rkdg_swe_kernels_processor.hpp"

namespace SWE {
namespace RKDG {
template <typename OMPISimUnitType>
void Problem::ompi_stage_kernel(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

    n_threads = (uint)omp_get_num_threads();
    thread_id = (uint)omp_get_thread_num();

    sim_per_thread = (sim_units.size() + n_threads - 1) / n_threads;

    begin_sim_id = sim_per_thread * thread_id;
    end_sim_id   = std::min(sim_per_thread * (thread_id + 1), (uint)sim_units.size());

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Current (time, stage): (" << sim_units[su_id]->stepper.GetTimeAtCurrentStage() << ','
                << sim_units[su_id]->stepper.GetStage() << ')' << std::endl;

            sim_units[su_id]->writer.GetLogFile() << "Exchanging data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(SWE::CommTypes::processor, sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::distributed_boundary_send_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->communicator.SendAll(SWE::CommTypes::processor, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::volume_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::source_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [&sim_units, su_id](auto& intface) { Problem::interface_kernel(sim_units[su_id]->stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&sim_units, su_id](auto& bound) { Problem::boundary_kernel(sim_units[su_id]->stepper, bound); });

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

        sim_units[su_id]->communicator.WaitAllReceives(SWE::CommTypes::processor,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::update_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = Problem::scrutinize_solution_kernel(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->communicator.WaitAllSends(SWE::CommTypes::processor,
                                                    sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Exchanging postprocessor data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(SWE::CommTypes::postprocessor,
                                                  sim_units[su_id]->stepper.GetTimestamp());

        if (SWE::PostProcessing::slope_limiting) {
            if (SWE::PostProcessing::wetting_drying) {
                sim_units[su_id]->discretization.mesh.CallForEachElement(
                    [&sim_units, su_id](auto& elt) { Problem::wetting_drying_kernel(sim_units[su_id]->stepper, elt); });
            }

            sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
                Problem::slope_limiting_prepare_element_kernel(sim_units[su_id]->stepper, elt);
            });

            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
                Problem::slope_limiting_distributed_boundary_send_kernel(sim_units[su_id]->stepper, dbound);
            });
        }

        sim_units[su_id]->communicator.SendAll(SWE::CommTypes::postprocessor, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting postprocessor work before receive" << std::endl;
        }

        if (SWE::PostProcessing::slope_limiting) {
            sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
                Problem::slope_limiting_prepare_interface_kernel(sim_units[su_id]->stepper, intface);
            });

            sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
                Problem::slope_limiting_prepare_boundary_kernel(sim_units[su_id]->stepper, bound);
            });
        }

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished postprocessor work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting to wait on postprocessor receive with timestamp: "
                                                  << sim_units[su_id]->stepper.GetTimestamp() << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(SWE::CommTypes::postprocessor,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting postprocessor work after receive" << std::endl;
        }

        if (SWE::PostProcessing::slope_limiting) {
            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
                Problem::slope_limiting_prepare_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
            });

            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::slope_limiting_kernel(sim_units[su_id]->stepper, elt); });
        }

        if (SWE::PostProcessing::wetting_drying) {
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::wetting_drying_kernel(sim_units[su_id]->stepper, elt); });
        }

        ++(sim_units[su_id]->stepper);

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished postprocessor work after receive" << std::endl
                                                  << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->communicator.WaitAllSends(SWE::CommTypes::postprocessor,
                                                    sim_units[su_id]->stepper.GetTimestamp());
    }
}
}
}

#endif