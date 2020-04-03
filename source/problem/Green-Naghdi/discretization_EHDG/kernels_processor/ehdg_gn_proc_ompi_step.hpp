#ifndef EHDG_GN_PROC_OMPI_STEP_HPP
#define EHDG_GN_PROC_OMPI_STEP_HPP

#include "ehdg_gn_kernels_processor.hpp"
#include "ehdg_gn_proc_ompi_sol_glob_prob.hpp"
#include "ehdg_gn_proc_derivatives_ompi.hpp"

namespace GN {
namespace EHDG {
template <template <typename> class OMPISimUnitType, typename ProblemType>
void seabed_update(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                   typename ProblemType::ProblemGlobalDataType& global_data,
                   typename ProblemType::ProblemStepperType& stepper,
                   const uint begin_sim_id,
                   const uint end_sim_id);

template <typename OMPISimUnitType>
void Problem::step_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                        ProblemGlobalDataType& global_data,
                        ProblemStepperType& stepper,
                        const uint begin_sim_id,
                        const uint end_sim_id) {
    auto& first_stepper  = stepper.GetFirstStepper();
    auto& second_stepper = stepper.GetSecondStepper();

    for (uint stage = 0; stage < first_stepper.GetNumStages(); ++stage) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (sim_units[su_id]->parser.ParsingInput()) {
                sim_units[su_id]->parser.ParseInput(first_stepper, sim_units[su_id]->discretization.mesh);
            }
        }
#pragma omp master
        { PetscLogStagePush(global_data.swe_stage); }

        SWE_SIM::Problem::stage_ompi(sim_units, global_data, first_stepper, begin_sim_id, end_sim_id);

#pragma omp master
        { PetscLogStagePop(); }
    }

    for (uint stage = 0; stage < second_stepper.GetNumStages(); ++stage) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (sim_units[su_id]->parser.ParsingInput()) {
                sim_units[su_id]->parser.ParseInput(second_stepper, sim_units[su_id]->discretization.mesh);
            }
        }

        Problem::dispersive_correction_ompi(sim_units, global_data, second_stepper, begin_sim_id, end_sim_id);
    }

    for (uint stage = 0; stage < first_stepper.GetNumStages(); ++stage) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            if (sim_units[su_id]->parser.ParsingInput()) {
                sim_units[su_id]->parser.ParseInput(first_stepper, sim_units[su_id]->discretization.mesh);
            }
        }
#pragma omp master
        { PetscLogStagePush(global_data.swe_stage); }

        SWE_SIM::Problem::stage_ompi(sim_units, global_data, first_stepper, begin_sim_id, end_sim_id);

#pragma omp master
        { PetscLogStagePop(); }
    }

    if (SWE::PostProcessing::bed_update && (second_stepper.GetStep() % SWE::PostProcessing::bed_update_freq == 0)) {
        seabed_update(sim_units, global_data, stepper, begin_sim_id, end_sim_id);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingOutput()) {
            sim_units[su_id]->writer.WriteOutput(stepper, sim_units[su_id]->discretization.mesh);
        }
    }
}

template <typename OMPISimUnitType>
void Problem::dispersive_correction_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                         ProblemGlobalDataType& global_data,
                                         ESSPRKStepper& stepper,
                                         const uint begin_sim_id,
                                         const uint end_sim_id) {
#pragma omp master
    { PetscLogStagePush(global_data.d_stage); }

    Problem::compute_derivatives_ompi(sim_units, global_data, stepper, begin_sim_id, end_sim_id);

#pragma omp master
    { PetscLogStagePop(); }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::local_dc_volume_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { Problem::local_dc_source_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&stepper](auto& edge_int) { Problem::local_dc_edge_interface_kernel(stepper, edge_int); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&stepper](auto& edge_bound) { Problem::local_dc_edge_boundary_kernel(stepper, edge_bound); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&stepper](auto& edge_dbound) { Problem::local_dc_edge_distributed_kernel(stepper, edge_dbound); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&stepper](auto& edge_int) { Problem::global_dc_edge_interface_kernel(stepper, edge_int); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&stepper](auto& edge_bound) { Problem::global_dc_edge_boundary_kernel(stepper, edge_bound); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&stepper](auto& edge_dbound) { Problem::global_dc_edge_distributed_kernel(stepper, edge_dbound); });
    }

    Problem::ompi_solve_global_dc_problem(sim_units, global_data, stepper, begin_sim_id, end_sim_id);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            Problem::dispersive_correction_kernel(stepper, elt);

            auto& state = elt.data.state[stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            stepper.UpdateState(elt);
        });
    }

#pragma omp barrier
#pragma omp master
    { ++(stepper); }
#pragma omp barrier
}

template <template <typename> class OMPISimUnitType, typename ProblemType>
void seabed_update(std::vector<std::unique_ptr<OMPISimUnitType<ProblemType>>>& sim_units,
                   typename ProblemType::ProblemGlobalDataType& global_data,
                   typename ProblemType::ProblemStepperType& stepper,
                   const uint begin_sim_id,
                   const uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(SWE_SIM::CommTypes::bound_state,
                                                  stepper.GetSecondStepper().GetTimestamp());
        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            SWE_SIM::Problem::distributed_boundary_send_kernel(stepper.GetFirstStepper(), dbound);
        });
        sim_units[su_id]->communicator.SendAll(SWE_SIM::CommTypes::bound_state,
                                               stepper.GetSecondStepper().GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(SWE_SIM::CommTypes::bound_state,
                                                       stepper.GetSecondStepper().GetTimestamp());
        SWE::seabed_update(stepper.GetSecondStepper(), sim_units[su_id]->discretization);
    }

    if (SWE::PostProcessing::bed_slope_limiting) {
        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->communicator.ReceiveAll(SWE_SIM::CommTypes::baryctr_state,
                                                      stepper.GetSecondStepper().GetTimestamp());
            sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
                auto& state              = elt.data.state[0];
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
                dbound.boundary_condition.exchanger.SetToSendBuffer(SWE_SIM::CommTypes::baryctr_state, message);
            });
            sim_units[su_id]->communicator.SendAll(SWE_SIM::CommTypes::baryctr_state,
                                                   stepper.GetSecondStepper().GetTimestamp());
        }

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->communicator.WaitAllReceives(SWE_SIM::CommTypes::baryctr_state,
                                                           stepper.GetSecondStepper().GetTimestamp());
            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
                std::vector<double> message(2);
                dbound.boundary_condition.exchanger.GetFromReceiveBuffer(SWE_SIM::CommTypes::baryctr_state, message);
                dbound.data.slope_limit_state.wet_neigh[dbound.bound_id]             = (bool)message[0];
                dbound.data.slope_limit_state.bath_at_baryctr_neigh[dbound.bound_id] = message[1];
            });
            SWE::CS_seabed_slope_limiter(sim_units[su_id]->discretization);
        }

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
            sim_units[su_id]->communicator.WaitAllSends(SWE_SIM::CommTypes::baryctr_state,
                                                        stepper.GetSecondStepper().GetTimestamp());
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        SWE::seabed_data_update(sim_units[su_id]->discretization);
        sim_units[su_id]->communicator.WaitAllSends(SWE_SIM::CommTypes::bound_state,
                                                    stepper.GetSecondStepper().GetTimestamp());
    }
    
    Problem::compute_bathymetry_derivatives_ompi(sim_units, global_data, begin_sim_id, end_sim_id);
}
}
}

#endif
