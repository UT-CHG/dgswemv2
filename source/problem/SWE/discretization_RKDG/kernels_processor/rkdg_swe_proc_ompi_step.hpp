#ifndef RKDG_SWE_PROC_OMPI_STEP_HPP
#define RKDG_SWE_PROC_OMPI_STEP_HPP

#include "rkdg_swe_kernels_processor.hpp"
#include "problem/SWE/problem_slope_limiter/swe_CS_sl_ompi.hpp"

namespace SWE {
namespace RKDG {
template <typename OMPISimType>
void Problem::step_ompi(OMPISimType* sim, uint begin_sim_id, uint end_sim_id) {
    assert(sim->sim_units.size() > 0);
    for (uint stage = 0; stage < sim->sim_units[0]->stepper.GetNumStages(); ++stage) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (sim->sim_units[su_id]->parser.ParsingInput()) {
                sim->sim_units[su_id]->parser.ParseInput(sim->sim_units[su_id]->stepper,
                                                         sim->sim_units[su_id]->discretization.mesh);
            }
        }

        Problem::stage_ompi(sim->sim_units, begin_sim_id, end_sim_id);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim->sim_units[su_id]->writer.WritingOutput()) {
            sim->sim_units[su_id]->writer.WriteOutput(sim->sim_units[su_id]->stepper,
                                                      sim->sim_units[su_id]->discretization.mesh);
        }
    }
}

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
            Problem::distributed_boundary_send_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::bound_state, sim_units[su_id]->stepper.GetTimestamp());
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

        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::bound_state,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            Problem::distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            auto& state = elt.data.state[sim_units[su_id]->stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            sim_units[su_id]->stepper.UpdateState(elt);
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });

        if (SWE::PostProcessing::wetting_drying) {
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::wetting_drying_kernel(sim_units[su_id]->stepper, elt); });
        }

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::bound_state, sim_units[su_id]->stepper.GetTimestamp());
    }

    if (SWE::PostProcessing::slope_limiting) {
        CS_slope_limiter_ompi(sim_units, begin_sim_id, end_sim_id, CommTypes::baryctr_state);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        ++(sim_units[su_id]->stepper);
    }
}
}
}

#endif