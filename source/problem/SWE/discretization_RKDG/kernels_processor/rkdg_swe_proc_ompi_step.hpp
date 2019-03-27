#ifndef RKDG_SWE_PROC_OMPI_STEP_HPP
#define RKDG_SWE_PROC_OMPI_STEP_HPP

#include "rkdg_swe_kernels_processor.hpp"
#include "problem/SWE/problem_slope_limiter/swe_CS_sl_ompi.hpp"

namespace SWE {
namespace RKDG {
template <template <typename> typename OMPISimUnitType, typename ProblemType>
void Problem::step_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                        typename ProblemType::ProblemGlobalDataType& global_data,
                        typename ProblemType::ProblemStepperType& stepper,
                        const uint begin_sim_id,
                        const uint end_sim_id) {
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
                         const uint begin_sim_id,
                         const uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Current (time, stage): (" << stepper.GetTimeAtCurrentStage()
                                                  << ',' << stepper.GetStage() << ')' << std::endl;

            sim_units[su_id]->writer.GetLogFile() << "Exchanging data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(CommTypes::bound_state, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [&stepper](auto& dbound) { Problem::distributed_boundary_send_kernel(stepper, dbound); });

        sim_units[su_id]->communicator.SendAll(CommTypes::bound_state, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::volume_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::source_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { Problem::interface_kernel(stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { Problem::boundary_kernel(stepper, bound); });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Starting to wait on receive with timestamp: " << stepper.GetTimestamp() << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::bound_state, stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [&stepper](auto& dbound) { Problem::distributed_boundary_kernel(stepper, dbound); });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            auto& state = elt.data.state[stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            stepper.UpdateState(elt);
        });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

#pragma omp barrier
#pragma omp master
    { ++(stepper); }
#pragma omp barrier

    if (SWE::PostProcessing::wetting_drying) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&stepper](auto& elt) { Problem::wetting_drying_kernel(stepper, elt); });
        }
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

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::bound_state, stepper.GetTimestamp());
    }
}
}
}

#endif