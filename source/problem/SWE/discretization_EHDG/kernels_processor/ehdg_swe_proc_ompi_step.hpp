#ifndef EHDG_SWE_PROC_OMPI_STEP_HPP
#define EHDG_SWE_PROC_OMPI_STEP_HPP

#include "ehdg_swe_kernels_processor.hpp"

namespace SWE {
namespace EHDG {
template <template <typename> class OMPISimUnitType, typename ProblemType>
void Problem::step_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                        typename ProblemType::ProblemGlobalDataType& global_data,
                        ProblemStepperType& stepper,
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

template <template <typename> class OMPISimUnitType, typename ProblemType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                         typename ProblemType::ProblemGlobalDataType& global_data,
                         ProblemStepperType& stepper,
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
            [&stepper](auto& dbound) { Problem::global_distributed_boundary_kernel(stepper, dbound); });

        sim_units[su_id]->communicator.SendAll(CommTypes::bound_state, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        /* Global Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { Problem::global_interface_kernel(stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { Problem::global_boundary_kernel(stepper, bound); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [](auto& edge_int) { edge_int.interface.specialization.ComputeNumericalFlux(edge_int); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&stepper](auto& edge_bound) {
            edge_bound.boundary.boundary_condition.ComputeNumericalFlux(stepper, edge_bound);
        });

        /* Global Pre Receive Step */

        /* Local Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::local_volume_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::local_source_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { Problem::local_interface_kernel(stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { Problem::local_boundary_kernel(stepper, bound); });
        /* Local Pre Receive Step */

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

        /* Global Post Receive Step */
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [](auto& edge_dbound) { edge_dbound.boundary.boundary_condition.ComputeNumericalFlux(edge_dbound); });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [&stepper](auto& dbound) { Problem::local_distributed_boundary_kernel(stepper, dbound); });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            auto& state = elt.data.state[stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            stepper.UpdateState(elt);
        });
        /* Local Post Receive Step */

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

    if (SWE::SedimentTransport::bed_update) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            SWE::seabed_update(stepper, sim_units[su_id]->discretization);
        }
    }

#pragma omp barrier
#pragma omp master
    { ++(stepper); }
#pragma omp barrier

    if (SWE::SedimentTransport::bed_slope_limiting) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->communicator.ReceiveAll(CommTypes::baryctr_state, stepper.GetTimestamp());
            sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
                auto& state              = elt.data.state[stepper.GetStage()];
                auto& sl_state           = elt.data.slope_limit_state;
                sl_state.bath_lin        = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
                sl_state.bath_at_baryctr = elt.ComputeLinearUbaryctr(sl_state.bath_lin);
            });
            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
                std::vector<double> message(2);
                message[0] = (double)dbound.data.wet_dry_state.wet;
                if (dbound.data.wet_dry_state.wet) {
                    message[1] = dbound.data.slope_limit_state.bath_at_baryctr;
                }
                dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::baryctr_state, message);
            });
            sim_units[su_id]->communicator.SendAll(CommTypes::baryctr_state, stepper.GetTimestamp());
        }

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->communicator.WaitAllReceives(CommTypes::baryctr_state, stepper.GetTimestamp());
            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
                std::vector<double> message(2);
                dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::baryctr_state, message);
                dbound.data.slope_limit_state.wet_neigh[dbound.bound_id]             = (bool)message[0];
                dbound.data.slope_limit_state.bath_at_baryctr_neigh[dbound.bound_id] = message[1];
            });
            SWE::CS_seabed_slope_limiter(stepper, sim_units[su_id]->discretization);
        }

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->communicator.WaitAllSends(CommTypes::baryctr_state, stepper.GetTimestamp());
        }
    }

    if (SWE::SedimentTransport::bed_update) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            SWE::seabed_data_update(stepper, sim_units[su_id]->discretization);
        }
    }

    if (SWE::PostProcessing::wetting_drying) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&stepper](auto& elt) { wetting_drying_kernel(stepper, elt); });
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
