#ifndef EHDG_GN_PROC_OMPI_STEP_HPP
#define EHDG_GN_PROC_OMPI_STEP_HPP

#include "general_definitions.hpp"

#include "ehdg_gn_kernels_processor.hpp"
#include "ehdg_gn_proc_ompi_sol_glob_prob.hpp"
#include "ehdg_gn_proc_derivatives_ompi.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimType>
void Problem::step_ompi(OMPISimType* sim, uint begin_sim_id, uint end_sim_id) {
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
        sim->sim_units[su_id]->discretization.mesh.CallForEachElement([sim, su_id](auto& elt) {
            uint n_stages = sim->sim_units[su_id]->stepper.GetNumStages();

            auto& state = elt.data.state;

            std::swap(state[0].q, state[n_stages].q);
        });

        if (sim->sim_units[su_id]->writer.WritingOutput()) {
            sim->sim_units[su_id]->writer.WriteOutput(sim->sim_units[su_id]->stepper,
                                                      sim->sim_units[su_id]->discretization.mesh);
        }
    }
}

template <typename OMPISimUnitType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, uint begin_sim_id, uint end_sim_id) {
    Problem::swe_stage_ompi(sim_units, begin_sim_id, end_sim_id);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& curr_state = elt.data.state[stage];
            auto& next_state = elt.data.state[stage + 1];

            curr_state.q = next_state.q;
        });
    }

    Problem::dispersive_correction_ompi(sim_units, begin_sim_id, end_sim_id);

    Problem::swe_stage_ompi(sim_units, begin_sim_id, end_sim_id);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });

        ++(sim_units[su_id]->stepper);
    }
}

template <typename OMPISimUnitType>
void Problem::swe_stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                             uint begin_sim_id,
                             uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Current (time, stage): (" << sim_units[su_id]->stepper.GetTimeAtCurrentStage() << ','
                << sim_units[su_id]->stepper.GetStage() << ')' << std::endl;

            sim_units[su_id]->writer.GetLogFile() << "Exchanging data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(SWE::EHDG::CommTypes::bound_state,
                                                  sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            SWE::EHDG::Problem::global_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->communicator.SendAll(SWE::EHDG::CommTypes::bound_state,
                                               sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work before receive" << std::endl;
        }

        /* Global Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            SWE::EHDG::Problem::global_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
            SWE::EHDG::Problem::global_boundary_kernel(sim_units[su_id]->stepper, bound);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&sim_units, su_id](auto& edge_int) {
            SWE::EHDG::Problem::global_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&sim_units, su_id](auto& edge_bound) {
            SWE::EHDG::Problem::global_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
        });
        /* Global Pre Receive Step */

        /* Local Pre Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            SWE::EHDG::Problem::local_volume_kernel(sim_units[su_id]->stepper, elt);
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            SWE::EHDG::Problem::local_source_kernel(sim_units[su_id]->stepper, elt);
        });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            SWE::EHDG::Problem::local_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
            SWE::EHDG::Problem::local_boundary_kernel(sim_units[su_id]->stepper, bound);
        });
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

        sim_units[su_id]->communicator.WaitAllReceives(SWE::EHDG::CommTypes::bound_state,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        /* Global Post Receive Step */
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&sim_units, su_id](auto& edge_dbound) {
                SWE::EHDG::Problem::global_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
            });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            SWE::EHDG::Problem::local_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& curr_state = elt.data.state[stage];
            auto& next_state = elt.data.state[stage + 1];

            curr_state.solution = elt.ApplyMinv(curr_state.rhs);

            next_state.q = curr_state.q + sim_units[su_id]->stepper.GetDT() / 2.0 * curr_state.solution;
        });
        /* Local Post Receive Step */

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->communicator.WaitAllSends(SWE::EHDG::CommTypes::bound_state,
                                                    sim_units[su_id]->stepper.GetTimestamp());
    }
}

template <typename OMPISimUnitType>
void Problem::dispersive_correction_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                         uint begin_sim_id,
                                         uint end_sim_id) {
    Problem::compute_derivatives_ompi(sim_units, begin_sim_id, end_sim_id);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::local_dc_volume_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { Problem::local_dc_source_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&sim_units, su_id](auto& edge_int) {
            Problem::local_dc_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&sim_units, su_id](auto& edge_bound) {
            Problem::local_dc_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&sim_units, su_id](auto& edge_dbound) {
                Problem::local_dc_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&sim_units, su_id](auto& edge_int) {
            Problem::global_dc_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&sim_units, su_id](auto& edge_bound) {
            Problem::global_dc_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&sim_units, su_id](auto& edge_dbound) {
                Problem::global_dc_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
            });
    }

    Problem::ompi_solve_global_dc_problem(sim_units, begin_sim_id, end_sim_id);

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            Problem::dispersive_correction_kernel(sim_units[su_id]->stepper, elt);

            auto& state = elt.data.state[sim_units[su_id]->stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            sim_units[su_id]->stepper.UpdateState(elt);
        });
    }
}
}
}

#endif