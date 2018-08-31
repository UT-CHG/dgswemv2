#ifndef IHDG_SWE_PROC_OMPI_STAGE_HPP
#define IHDG_SWE_PROC_OMPI_STAGE_HPP

#include "general_definitions.hpp"

#include "ihdg_swe_kernels_processor.hpp"

namespace SWE {
namespace IHDG {
template <typename OMPISimUnitType>
void Problem::ompi_stage_kernel(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

    n_threads = (uint)omp_get_num_threads();
    thread_id = (uint)omp_get_thread_num();

    sim_per_thread = (sim_units.size() + n_threads - 1) / n_threads;

    begin_sim_id = sim_per_thread * thread_id;
    end_sim_id   = std::min(sim_per_thread * (thread_id + 1), (uint)sim_units.size());

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
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

        for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
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

        if (iter == 10) {
            break;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = Problem::scrutinize_solution_kernel(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });

        ++(sim_units[su_id]->stepper);
    }
}
}
}

#endif