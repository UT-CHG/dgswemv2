#ifndef EHDG_GN_PROC_OMPI_SOL_GLOB_PROB_HPP
#define EHDG_GN_PROC_OMPI_SOL_GLOB_PROB_HPP

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::ompi_solve_global_dc_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                           ProblemGlobalDataType& global_data,
                                           const ESSPRKStepper& stepper,
                                           const uint begin_sim_id,
                                           const uint end_sim_id) {
#pragma omp master
    { PetscLogStagePush(global_data.con_stage); }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            if (elt.data.wet_dry_state.wet /*&& elt.data.source.dispersive_correction*/) {
                auto& internal     = elt.data.internal;
                internal.w2_w2_inv = inverse(internal.w2_w2);
                internal.w1_w1 -= internal.w1_w2 * internal.w2_w2_inv * internal.w2_w1;
                for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
                    elt.data.boundary[bound_id].w1_w1_hat -=
                        internal.w1_w2 * internal.w2_w2_inv * elt.data.boundary[bound_id].w2_w1_hat;
                    solve_sle(internal.w1_w1, elt.data.boundary[bound_id].w1_w1_hat);
                }
                solve_sle(internal.w1_w1, internal.w1_rhs);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([](auto& edge_int) {
            if ((edge_int.interface.data_in.wet_dry_state.wet /*&&
                 edge_int.interface.data_in.source.dispersive_correction*/) ||
                (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                 edge_int.interface.data_ex.source.dispersive_correction*/)) {
                auto& edge_internal = edge_int.edge_data.edge_internal;
                auto& internal_in   = edge_int.interface.data_in.internal;
                auto& internal_ex   = edge_int.interface.data_ex.internal;
                auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
                auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

                set_constant(edge_internal.w1_hat_rhs, 0.0);

                if (edge_int.interface.data_in.wet_dry_state.wet /*&&
                    edge_int.interface.data_in.source.dispersive_correction*/) {
                    boundary_in.w1_hat_w1 -= boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * internal_in.w2_w1;
                    edge_internal.w1_hat_w1_hat -=
                        boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_in.w2_w1_hat +
                        boundary_in.w1_hat_w1 * boundary_in.w1_w1_hat;
                    edge_internal.w1_hat_rhs -= boundary_in.w1_hat_w1 * internal_in.w1_rhs;
                }

                if (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                    edge_int.interface.data_ex.source.dispersive_correction*/) {
                    boundary_ex.w1_hat_w1 -= boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * internal_ex.w2_w1;
                    edge_internal.w1_hat_w1_hat -=
                        boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * boundary_ex.w2_w1_hat +
                        boundary_ex.w1_hat_w1 * boundary_ex.w1_w1_hat;
                    edge_internal.w1_hat_rhs -= boundary_ex.w1_hat_w1 * internal_ex.w1_rhs;
                }

                edge_internal.w1_hat_w1_hat_flat = flatten<double>(edge_internal.w1_hat_w1_hat);

                uint bcon_id = 0;
                if (edge_int.interface.data_in.wet_dry_state.wet /*&&
                    edge_int.interface.data_in.source.dispersive_correction*/) {
                    for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
                        if (bound_id == edge_int.interface.bound_id_in)
                            continue;
                        auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];
                        edge_internal.w1_hat_w1_hat =
                            -(boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_con.w2_w1_hat +
                              boundary_in.w1_hat_w1 * boundary_con.w1_w1_hat);
                        edge_internal.w1_hat_w1_hat_con_flat[bcon_id] = flatten<double>(edge_internal.w1_hat_w1_hat);
                        ++bcon_id;
                    }
                }

                if (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                    edge_int.interface.data_ex.source.dispersive_correction*/) {
                    for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
                        if (bound_id == edge_int.interface.bound_id_ex)
                            continue;
                        auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];
                        edge_internal.w1_hat_w1_hat =
                            -(boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * boundary_con.w2_w1_hat +
                              boundary_ex.w1_hat_w1 * boundary_con.w1_w1_hat);
                        edge_internal.w1_hat_w1_hat_con_flat[bcon_id] = flatten<double>(edge_internal.w1_hat_w1_hat);
                        ++bcon_id;
                    }
                }
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([](auto& edge_bound) {
            if (edge_bound.boundary.data.wet_dry_state.wet /*&&
                edge_bound.boundary.data.source.dispersive_correction*/) {
                auto& edge_internal = edge_bound.edge_data.edge_internal;
                auto& internal      = edge_bound.boundary.data.internal;
                auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

                boundary.w1_hat_w1 -= boundary.w1_hat_w2 * internal.w2_w2_inv * internal.w2_w1;
                edge_internal.w1_hat_w1_hat -= boundary.w1_hat_w2 * internal.w2_w2_inv * boundary.w2_w1_hat +
                                               boundary.w1_hat_w1 * boundary.w1_w1_hat;
                edge_internal.w1_hat_rhs         = -boundary.w1_hat_w1 * internal.w1_rhs;
                edge_internal.w1_hat_w1_hat_flat = flatten<double>(edge_internal.w1_hat_w1_hat);

                uint bcon_id = 0;
                for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
                    if (bound_id == edge_bound.boundary.bound_id)
                        continue;
                    auto& boundary_con          = edge_bound.boundary.data.boundary[bound_id];
                    edge_internal.w1_hat_w1_hat = -(boundary.w1_hat_w2 * internal.w2_w2_inv * boundary_con.w2_w1_hat +
                                                    boundary.w1_hat_w1 * boundary_con.w1_w1_hat);
                    edge_internal.w1_hat_w1_hat_con_flat[bcon_id] = flatten<double>(edge_internal.w1_hat_w1_hat);
                    ++bcon_id;
                }
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([](auto& edge_dbound) {
            if (edge_dbound.boundary.data.wet_dry_state.wet /*&&
                edge_dbound.boundary.data.source.dispersive_correction*/) {
                auto& edge_internal = edge_dbound.edge_data.edge_internal;
                auto& internal      = edge_dbound.boundary.data.internal;
                auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

                boundary.w1_hat_w1 -= boundary.w1_hat_w2 * internal.w2_w2_inv * internal.w2_w1;
                edge_internal.w1_hat_w1_hat -= boundary.w1_hat_w2 * internal.w2_w2_inv * boundary.w2_w1_hat +
                                               boundary.w1_hat_w1 * boundary.w1_w1_hat;
                edge_internal.w1_hat_rhs         = -boundary.w1_hat_w1 * internal.w1_rhs;
                edge_internal.w1_hat_w1_hat_flat = flatten<double>(edge_internal.w1_hat_w1_hat);

                uint bcon_id = 0;
                for (uint bound_id = 0; bound_id < edge_dbound.boundary.data.get_nbound(); ++bound_id) {
                    if (bound_id == edge_dbound.boundary.bound_id)
                        continue;
                    auto& boundary_con          = edge_dbound.boundary.data.boundary[bound_id];
                    edge_internal.w1_hat_w1_hat = -(boundary.w1_hat_w2 * internal.w2_w2_inv * boundary_con.w2_w1_hat +
                                                    boundary.w1_hat_w1 * boundary_con.w1_w1_hat);
                    edge_internal.w1_hat_w1_hat_con_flat[bcon_id] = flatten<double>(edge_internal.w1_hat_w1_hat);
                    ++bcon_id;
                }
            }
        });
    }

#pragma omp barrier
#pragma omp master
    {
        if (SWE::PostProcessing::wetting_drying)
            reset_PETSC_solver(sim_units, global_data, stepper, begin_sim_id, end_sim_id);

        Mat& w1_hat_w1_hat = global_data.w1_hat_w1_hat;
        Vec& w1_hat_rhs    = global_data.w1_hat_rhs;

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&w1_hat_rhs, &w1_hat_w1_hat](auto& edge_int) {
                    if ((edge_int.interface.data_in.wet_dry_state.wet /*&&
                         edge_int.interface.data_in.source.dispersive_correction*/) ||
                            (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                         edge_int.interface.data_ex.source.dispersive_correction*/)) {
                        auto& edge_internal = edge_int.edge_data.edge_internal;

                        VecSetValuesBlocked(w1_hat_rhs,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            edge_internal.w1_hat_rhs.data(),
                                            ADD_VALUES);
                        MatSetValuesBlocked(w1_hat_w1_hat,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            edge_internal.w1_hat_w1_hat_flat.data(),
                                            ADD_VALUES);

                        uint bcon_id = 0;
                        if (edge_int.interface.data_in.wet_dry_state.wet /*&&
                            edge_int.interface.data_in.source.dispersive_correction*/) {
                            for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
                                if (bound_id == edge_int.interface.bound_id_in)
                                    continue;
                                auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];
                                MatSetValuesBlocked(w1_hat_w1_hat,
                                                    1,
                                                    (int*)&edge_internal.dc_global_dof_indx,
                                                    1,
                                                    (int*)&boundary_con.dc_global_dof_indx,
                                                    edge_internal.w1_hat_w1_hat_con_flat[bcon_id].data(),
                                                    ADD_VALUES);
                                ++bcon_id;
                            }
                        }

                        if (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                            edge_int.interface.data_ex.source.dispersive_correction*/) {
                            for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
                                if (bound_id == edge_int.interface.bound_id_ex)
                                    continue;
                                auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];
                                MatSetValuesBlocked(w1_hat_w1_hat,
                                                    1,
                                                    (int*)&edge_internal.dc_global_dof_indx,
                                                    1,
                                                    (int*)&boundary_con.dc_global_dof_indx,
                                                    edge_internal.w1_hat_w1_hat_con_flat[bcon_id].data(),
                                                    ADD_VALUES);
                                ++bcon_id;
                            }
                        }
                    }
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&w1_hat_rhs, &w1_hat_w1_hat](auto& edge_bound) {
                    if (edge_bound.boundary.data.wet_dry_state.wet /*&&
                        edge_bound.boundary.data.source.dispersive_correction*/) {
                        auto& edge_internal = edge_bound.edge_data.edge_internal;

                        VecSetValuesBlocked(w1_hat_rhs,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            edge_internal.w1_hat_rhs.data(),
                                            ADD_VALUES);
                        MatSetValuesBlocked(w1_hat_w1_hat,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            edge_internal.w1_hat_w1_hat_flat.data(),
                                            ADD_VALUES);

                        uint bcon_id = 0;
                        for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
                            if (bound_id == edge_bound.boundary.bound_id)
                                continue;
                            auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];
                            MatSetValuesBlocked(w1_hat_w1_hat,
                                                1,
                                                (int*)&edge_internal.dc_global_dof_indx,
                                                1,
                                                (int*)&boundary_con.dc_global_dof_indx,
                                                edge_internal.w1_hat_w1_hat_con_flat[bcon_id].data(),
                                                ADD_VALUES);
                            ++bcon_id;
                        }
                    }
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&w1_hat_rhs, &w1_hat_w1_hat](auto& edge_dbound) {
                    if (edge_dbound.boundary.data.wet_dry_state.wet /*&&
                        edge_dbound.boundary.data.source.dispersive_correction*/) {
                        auto& edge_internal = edge_dbound.edge_data.edge_internal;

                        VecSetValuesBlocked(w1_hat_rhs,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            edge_internal.w1_hat_rhs.data(),
                                            ADD_VALUES);
                        MatSetValuesBlocked(w1_hat_w1_hat,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            1,
                                            (int*)&edge_internal.dc_global_dof_indx,
                                            edge_internal.w1_hat_w1_hat_flat.data(),
                                            ADD_VALUES);

                        uint bcon_id = 0;
                        for (uint bound_id = 0; bound_id < edge_dbound.boundary.data.get_nbound(); ++bound_id) {
                            if (bound_id == edge_dbound.boundary.bound_id)
                                continue;
                            auto& boundary_con = edge_dbound.boundary.data.boundary[bound_id];
                            MatSetValuesBlocked(w1_hat_w1_hat,
                                                1,
                                                (int*)&edge_internal.dc_global_dof_indx,
                                                1,
                                                (int*)&boundary_con.dc_global_dof_indx,
                                                edge_internal.w1_hat_w1_hat_con_flat[bcon_id].data(),
                                                ADD_VALUES);
                            ++bcon_id;
                        }
                    }
                });
        }

        MatAssemblyBegin(w1_hat_w1_hat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(w1_hat_w1_hat, MAT_FINAL_ASSEMBLY);

        VecAssemblyBegin(w1_hat_rhs);
        VecAssemblyEnd(w1_hat_rhs);

        KSP& dc_ksp            = global_data.dc_ksp;
        Vec& dc_sol            = global_data.dc_sol;
        VecScatter& dc_scatter = global_data.dc_scatter;

        PetscLogStagePop();
        PetscLogStagePush(global_data.sol_stage);

        KSPSolve(dc_ksp, w1_hat_rhs, w1_hat_rhs);

        PetscLogStagePop();
        PetscLogStagePush(global_data.prop_stage);

        VecScatterBegin(dc_scatter, w1_hat_rhs, dc_sol, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(dc_scatter, w1_hat_rhs, dc_sol, INSERT_VALUES, SCATTER_FORWARD);

        double* sol_ptr;
        VecGetArray(dc_sol, &sol_ptr);

        int sol_size;
        VecGetLocalSize(dc_sol, &sol_size);

        global_data.dc_solution = vector_from_array(sol_ptr, sol_size);

        VecRestoreArray(dc_sol, &sol_ptr);

        MatZeroEntries(w1_hat_w1_hat);
        VecZeroEntries(w1_hat_rhs);

        if (stepper.GetStep() == 0 && !SWE::PostProcessing::wetting_drying) {
            MatSetOption(w1_hat_w1_hat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
        }
    }
#pragma omp barrier

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&global_data](auto& edge_int) {
            auto& edge_internal = edge_int.edge_data.edge_internal;
            auto& internal_in   = edge_int.interface.data_in.internal;
            auto& internal_ex   = edge_int.interface.data_ex.internal;
            auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
            auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

            if ((edge_int.interface.data_in.wet_dry_state.wet /*&&
                 edge_int.interface.data_in.source.dispersive_correction*/) ||
                (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                 edge_int.interface.data_ex.source.dispersive_correction*/)) {
                const uint n_global_dofs = edge_int.edge_data.get_ndof() * GN::n_dimensions;
                const auto w1_hat =
                    subvector(global_data.dc_solution, edge_internal.dc_local_dof_indx * n_global_dofs, n_global_dofs);
                if (edge_int.interface.data_in.wet_dry_state.wet /*&&
                edge_int.interface.data_in.source.dispersive_correction*/)
                    internal_in.w1_rhs -= boundary_in.w1_w1_hat * w1_hat;
                if (edge_int.interface.data_ex.wet_dry_state.wet /*&&
                edge_int.interface.data_ex.source.dispersive_correction*/)
                    internal_ex.w1_rhs -= boundary_ex.w1_w1_hat * w1_hat;
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_data](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;
            auto& internal      = edge_bound.boundary.data.internal;
            auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            if (edge_bound.boundary.data.wet_dry_state.wet /*&&
                edge_bound.boundary.data.source.dispersive_correction*/) {
                const uint n_global_dofs = edge_bound.edge_data.get_ndof() * GN::n_dimensions;
                const auto w1_hat =
                    subvector(global_data.dc_solution, edge_internal.dc_local_dof_indx * n_global_dofs, n_global_dofs);
                internal.w1_rhs -= boundary.w1_w1_hat * w1_hat;
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&global_data](auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;
            auto& internal      = edge_dbound.boundary.data.internal;
            auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            if (edge_dbound.boundary.data.wet_dry_state.wet /*&&
                edge_dbound.boundary.data.source.dispersive_correction*/) {
                const uint n_global_dofs = edge_dbound.edge_data.get_ndof() * GN::n_dimensions;
                const auto w1_hat =
                    subvector(global_data.dc_solution, edge_internal.dc_local_dof_indx * n_global_dofs, n_global_dofs);
                internal.w1_rhs -= boundary.w1_w1_hat * w1_hat;
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            const uint stage = stepper.GetStage();
            auto& state      = elt.data.state[stage];
            auto& internal   = elt.data.internal;
            if (elt.data.wet_dry_state.wet /*&& elt.data.source.dispersive_correction*/)
                state.w1 = reshape<double, GN::n_dimensions, SO::ColumnMajor>(internal.w1_rhs, elt.data.get_ndof());
        });
    }

#pragma omp master
    { PetscLogStagePop(); }
}

template <typename OMPISimUnitType>
void Problem::reset_PETSC_solver(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                 ProblemGlobalDataType& global_data,
                                 const ESSPRKStepper& stepper,
                                 const uint begin_sim_id,
                                 const uint end_sim_id) {
    MatDestroy(&global_data.w1_hat_w1_hat);
    VecDestroy(&global_data.w1_hat_rhs);
    KSPDestroy(&global_data.dc_ksp);

    ISDestroy(&global_data.dc_from);
    ISDestroy(&global_data.dc_to);
    VecScatterDestroy(&global_data.dc_scatter);
    VecDestroy(&global_data.dc_sol);

    uint dc_global_dof_offset = 0;
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&dc_global_dof_offset](auto& edge_int) {
                auto& edge_internal = edge_int.edge_data.edge_internal;
                auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
                auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

                // Set indexes for global matrix construction
                if ((edge_int.interface.data_in.wet_dry_state
                         .wet /*&& edge_int.interface.data_in.source.dispersive_correction*/) ||
                    (edge_int.interface.data_ex.wet_dry_state
                         .wet /*&& edge_int.interface.data_ex.source.dispersive_correction*/)) {
                    edge_internal.dc_global_dof_indx = dc_global_dof_offset;
                    boundary_in.dc_global_dof_indx   = edge_internal.dc_global_dof_indx;
                    boundary_ex.dc_global_dof_indx   = edge_internal.dc_global_dof_indx;
                    ++dc_global_dof_offset;
                }
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&dc_global_dof_offset](auto& edge_bound) {
                auto& edge_internal = edge_bound.edge_data.edge_internal;
                auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

                // Set indexes for global matrix construction
                if (edge_bound.boundary.data.wet_dry_state
                        .wet /*&& edge_bound.boundary.data.source.dispersive_correction*/) {
                    edge_internal.dc_global_dof_indx = dc_global_dof_offset;
                    boundary.dc_global_dof_indx      = edge_internal.dc_global_dof_indx;
                    ++dc_global_dof_offset;
                }
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&dc_global_dof_offset](
                                                                                      auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;
            auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            const uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
            const uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
            const uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
            const uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

            // Set indexes for global matrix construction
            if (edge_dbound.boundary.data.wet_dry_state.wet /*&& edge_bound.boundary.data.source.dispersive_correction*/
                || edge_dbound.boundary.boundary_condition.wet_neighbor) {
                if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
                    edge_internal.dc_global_dof_indx = dc_global_dof_offset;
                    boundary.dc_global_dof_indx      = edge_internal.dc_global_dof_indx;
                    ++dc_global_dof_offset;
                }
            }
        });
    }

    int n_localities;
    int locality_id;
    MPI_Comm_size(MPI_COMM_WORLD, &n_localities);
    MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

    std::vector<uint> dc_global_dof_offsets;
    if (locality_id == 0) {
        dc_global_dof_offsets.resize(n_localities);
    }
    MPI_Gather(
        &dc_global_dof_offset, 1, MPI_UNSIGNED, &dc_global_dof_offsets.front(), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    uint n_dc_global_dofs = 0;
    if (locality_id == 0) {
        n_dc_global_dofs = std::accumulate(dc_global_dof_offsets.begin(), dc_global_dof_offsets.end(), 0);
    }
    MPI_Bcast(&n_dc_global_dofs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    const uint n_global_dofs = (sim_units.front()->discretization.mesh.GetP() + 1) * GN::n_dimensions;
    MatCreateBAIJ(MPI_COMM_WORLD,
                  n_global_dofs,
                  n_global_dofs * dc_global_dof_offset,
                  n_global_dofs * dc_global_dof_offset,
                  n_global_dofs * n_dc_global_dofs,
                  n_global_dofs * n_dc_global_dofs,
                  5,
                  NULL,
                  4,
                  NULL,
                  &(global_data.w1_hat_w1_hat));

    VecCreateMPI(MPI_COMM_WORLD,
                 n_global_dofs * dc_global_dof_offset,
                 n_global_dofs * n_dc_global_dofs,
                 &(global_data.w1_hat_rhs));
    VecSetBlockSize(global_data.w1_hat_rhs, n_global_dofs);

    KSPCreate(MPI_COMM_WORLD, &(global_data.dc_ksp));
    KSPSetOperators(global_data.dc_ksp, global_data.w1_hat_w1_hat, global_data.w1_hat_w1_hat);
    // KSPGetPC(global_data.dc_ksp, &(global_data.dc_pc));
    // PCSetType(global_data.dc_pc, PCLU);

    if (locality_id == 0) {
        std::rotate(dc_global_dof_offsets.begin(), dc_global_dof_offsets.end() - 1, dc_global_dof_offsets.end());
        dc_global_dof_offsets.front() = 0;
        for (int locality_id = 1; locality_id < n_localities; ++locality_id) {
            dc_global_dof_offsets[locality_id] += dc_global_dof_offsets[locality_id - 1];
        }
    }
    MPI_Scatter(
        &dc_global_dof_offsets.front(), 1, MPI_UNSIGNED, &dc_global_dof_offset, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dc_global_dof_indx, stepper.GetTimestamp());
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&dc_global_dof_offset](auto& edge_int) {
                auto& edge_internal = edge_int.edge_data.edge_internal;
                auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
                auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

                if ((edge_int.interface.data_in.wet_dry_state
                         .wet /*&& edge_int.interface.data_in.source.dispersive_correction*/) ||
                    (edge_int.interface.data_ex.wet_dry_state
                         .wet /*&& edge_int.interface.data_ex.source.dispersive_correction*/)) {
                    edge_internal.dc_global_dof_indx += dc_global_dof_offset;
                    boundary_in.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
                    boundary_ex.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
                }
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&dc_global_dof_offset](auto& edge_bound) {
                auto& edge_internal = edge_bound.edge_data.edge_internal;
                auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

                if (edge_bound.boundary.data.wet_dry_state
                        .wet /*&& edge_bound.boundary.data.source.dispersive_correction*/) {
                    edge_internal.dc_global_dof_indx += dc_global_dof_offset;
                    boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
                }
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&dc_global_dof_offset](
                                                                                      auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;
            auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            const uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
            const uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
            const uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
            const uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

            if (edge_dbound.boundary.data.wet_dry_state.wet /*&& edge_bound.boundary.data.source.dispersive_correction*/
                || edge_dbound.boundary.boundary_condition.wet_neighbor) {
                if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
                    edge_internal.dc_global_dof_indx += dc_global_dof_offset;
                    boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
                    std::vector<double> message(1, (double)edge_internal.dc_global_dof_indx);
                    edge_dbound.boundary.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dc_global_dof_indx,
                                                                                      message);
                }
            }
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dc_global_dof_indx, stepper.GetTimestamp());
    }

    std::vector<uint> dc_global_dof_indx;
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dc_global_dof_indx, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&dc_global_dof_indx](auto& edge_int) {
            auto& edge_internal = edge_int.edge_data.edge_internal;

            if ((edge_int.interface.data_in.wet_dry_state
                     .wet /*&& edge_int.interface.data_in.source.dispersive_correction*/) ||
                (edge_int.interface.data_ex.wet_dry_state
                     .wet /*&& edge_int.interface.data_ex.source.dispersive_correction*/)) {
                edge_internal.dc_local_dof_indx = dc_global_dof_indx.size();
                dc_global_dof_indx.push_back(edge_internal.dc_global_dof_indx);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&dc_global_dof_indx](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            if (edge_bound.boundary.data.wet_dry_state
                    .wet /*&& edge_bound.boundary.data.source.dispersive_correction*/) {
                edge_internal.dc_local_dof_indx = dc_global_dof_indx.size();
                dc_global_dof_indx.push_back(edge_internal.dc_global_dof_indx);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&dc_global_dof_indx](auto& edge_dbound) {
                auto& edge_internal = edge_dbound.edge_data.edge_internal;
                auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

                const uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
                const uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
                const uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
                const uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

                if (edge_dbound.boundary.data.wet_dry_state
                        .wet /*&& edge_bound.boundary.data.source.dispersive_correction*/) {
                    if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
                        std::vector<double> message(1);
                        edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(
                            CommTypes::dc_global_dof_indx, message);
                        edge_internal.dc_global_dof_indx = (uint)message[0];
                        boundary.dc_global_dof_indx      = edge_internal.dc_global_dof_indx;
                    }

                    edge_internal.dc_local_dof_indx = dc_global_dof_indx.size();
                    dc_global_dof_indx.push_back(edge_internal.dc_global_dof_indx);
                }
            });
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dc_global_dof_indx, stepper.GetTimestamp());
    }

    VecCreateSeq(MPI_COMM_SELF, n_global_dofs * dc_global_dof_indx.size(), &(global_data.dc_sol));

    ISCreateBlock(MPI_COMM_SELF,
                  n_global_dofs,
                  dc_global_dof_indx.size(),
                  (int*)&dc_global_dof_indx.front(),
                  PETSC_COPY_VALUES,
                  &(global_data.dc_from));
    ISCreateStride(MPI_COMM_SELF, n_global_dofs * dc_global_dof_indx.size(), 0, 1, &(global_data.dc_to));

    VecScatterCreate(
        global_data.w1_hat_rhs, global_data.dc_from, global_data.dc_sol, global_data.dc_to, &(global_data.dc_scatter));
}
}
}

#endif
