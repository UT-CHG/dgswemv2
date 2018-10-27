#ifndef IHDG_SWE_PROC_OMPI_STEP_HPP
#define IHDG_SWE_PROC_OMPI_STEP_HPP

#include "ihdg_swe_kernels_processor.hpp"
#include "ihdg_swe_proc_ompi_sol_glob_prob.hpp"

namespace SWE {
namespace IHDG {
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
    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state      = elt.data.state[stage + 1];
            auto& state_prev = elt.data.state[stage];
            auto& internal   = elt.data.internal;

            state.q = state_prev.q;

            internal.q_prev_at_gp = elt.ComputeUgp(state_prev.q);
        });
    }

    uint iter = 0;
    while (true) {
        iter++;

        for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
            /* Local Step */
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::local_volume_kernel(sim_units[su_id]->stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::local_source_kernel(sim_units[su_id]->stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
                Problem::local_interface_kernel(sim_units[su_id]->stepper, intface);
            });

            sim_units[su_id]->discretization.mesh.CallForEachBoundary(
                [&sim_units, su_id](auto& bound) { Problem::local_boundary_kernel(sim_units[su_id]->stepper, bound); });

            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
                Problem::local_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
            });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&sim_units, su_id](auto& edge_int) {
                    Problem::local_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&sim_units, su_id](auto& edge_bound) {
                    Problem::local_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&sim_units, su_id](auto& edge_dbound) {
                    Problem::local_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
                });
            /* Local Step */

            /* Global Step */
            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&sim_units, su_id](auto& edge_int) {
                    Problem::global_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&sim_units, su_id](auto& edge_bound) {
                    Problem::global_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&sim_units, su_id](auto& edge_dbound) {
                    Problem::global_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
                });
            /* Global Step */
        }

        bool converged = ompi_solve_global_problem(sim_units, begin_sim_id, end_sim_id);

        if (converged) {
            break;
        }

        if (iter == 100) {
            break;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });

        ++(sim_units[su_id]->stepper);
    }
}
}
}

#endif