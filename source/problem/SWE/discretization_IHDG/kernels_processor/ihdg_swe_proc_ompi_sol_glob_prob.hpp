#ifndef IHDG_SWE_PROC_OMPI_SOL_GLOB_PROB_HPP
#define IHDG_SWE_PROC_OMPI_SOL_GLOB_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename OMPISimUnitType>
bool Problem::ompi_solve_global_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
#pragma omp master
    {
        uint local_dof_offset  = 3 * 2048;
        uint global_dof_offset = 2 * 3136;

        SparseMatrix<double> delta_local_inv;
        SparseMatrix<double> delta_hat_local;
        DynVector<double> rhs_local;

        SparseMatrix<double> delta_global;
        SparseMatrix<double> delta_hat_global;
        DynVector<double> rhs_global;

        delta_local_inv.resize(local_dof_offset * SWE::n_variables, local_dof_offset * SWE::n_variables);
        delta_hat_local.resize(local_dof_offset * SWE::n_variables, global_dof_offset * SWE::n_variables);
        rhs_local.resize(local_dof_offset * SWE::n_variables);

        delta_global.resize(global_dof_offset * SWE::n_variables, local_dof_offset * SWE::n_variables);
        delta_hat_global.resize(global_dof_offset * SWE::n_variables, global_dof_offset * SWE::n_variables);
        rhs_global.resize(global_dof_offset * SWE::n_variables);

        set_constant(rhs_local, 0.0);
        set_constant(rhs_global, 0.0);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            SparseMatrix<double> delta_local_inv_temp;
            SparseMatrix<double> delta_hat_local_temp;
            DynVector<double> rhs_local_temp;

            SparseMatrix<double> delta_global_temp;
            SparseMatrix<double> delta_hat_global_temp;
            DynVector<double> rhs_global_temp;

            delta_local_inv_temp.resize(local_dof_offset * SWE::n_variables, local_dof_offset * SWE::n_variables);
            delta_hat_local_temp.resize(local_dof_offset * SWE::n_variables, global_dof_offset * SWE::n_variables);
            rhs_local_temp.resize(local_dof_offset * SWE::n_variables);

            delta_global_temp.resize(global_dof_offset * SWE::n_variables, local_dof_offset * SWE::n_variables);
            delta_hat_global_temp.resize(global_dof_offset * SWE::n_variables, global_dof_offset * SWE::n_variables);
            rhs_global_temp.resize(global_dof_offset * SWE::n_variables);

            SparseMatrixMeta<double> sparse_delta_local_inv;
            SparseMatrixMeta<double> sparse_delta_hat_local;
            SparseMatrixMeta<double> sparse_delta_global;
            SparseMatrixMeta<double> sparse_delta_hat_global;

            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&rhs_local_temp, &sparse_delta_local_inv](auto& elt) {
                    auto& internal = elt.data.internal;

                    uint ndof         = elt.data.get_ndof();
                    uint local_offset = internal.local_dof_offset;

                    subvector(rhs_local_temp, local_offset * SWE::n_variables, ndof * SWE::n_variables) =
                        internal.rhs_local;

                    internal.delta_local_inv = inverse(internal.delta_local);

                    for (uint i = 0; i < ndof * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof * SWE::n_variables; ++j) {
                            sparse_delta_local_inv.add_triplet(local_offset * SWE::n_variables + i,
                                                               local_offset * SWE::n_variables + j,
                                                               internal.delta_local_inv(i, j));
                        }
                    }
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&rhs_global_temp, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](
                    auto& edge_int) {
                    auto& edge_internal = edge_int.edge_data.edge_internal;

                    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
                    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

                    uint ndof_local_in = edge_int.interface.data_in.get_ndof();
                    uint ndof_local_ex = edge_int.interface.data_ex.get_ndof();
                    uint ndof_global   = edge_int.edge_data.get_ndof();

                    uint local_offset_in = boundary_in.local_dof_offset;
                    uint local_offset_ex = boundary_ex.local_dof_offset;
                    uint global_offset   = edge_internal.global_dof_offset;

                    subvector(rhs_global_temp, global_offset * SWE::n_variables, ndof_global * SWE::n_variables) =
                        edge_internal.rhs_global;

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_local_in * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_local.add_triplet(local_offset_in * SWE::n_variables + i,
                                                               global_offset * SWE::n_variables + j,
                                                               boundary_in.delta_hat_local(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_local_ex * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_local.add_triplet(local_offset_ex * SWE::n_variables + i,
                                                               global_offset * SWE::n_variables + j,
                                                               boundary_ex.delta_hat_local(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_local_in * SWE::n_variables; ++j) {
                            sparse_delta_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            local_offset_in * SWE::n_variables + j,
                                                            boundary_in.delta_global(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_local_ex * SWE::n_variables; ++j) {
                            sparse_delta_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            local_offset_ex * SWE::n_variables + j,
                                                            boundary_ex.delta_global(i, j));
                        }
                    }
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&rhs_global_temp, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](
                    auto& edge_bound) {
                    auto& edge_internal = edge_bound.edge_data.edge_internal;

                    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

                    uint ndof_local  = edge_bound.boundary.data.get_ndof();
                    uint ndof_global = edge_bound.edge_data.get_ndof();

                    uint local_offset  = boundary.local_dof_offset;
                    uint global_offset = edge_internal.global_dof_offset;

                    subvector(rhs_global_temp, global_offset * SWE::n_variables, ndof_global * SWE::n_variables) =
                        edge_internal.rhs_global;

                    for (uint i = 0; i < ndof_local * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_local.add_triplet(local_offset * SWE::n_variables + i,
                                                               global_offset * SWE::n_variables + j,
                                                               boundary.delta_hat_local(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_local * SWE::n_variables; ++j) {
                            sparse_delta_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            local_offset * SWE::n_variables + j,
                                                            boundary.delta_global(i, j));
                        }
                    }
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&rhs_global_temp, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](
                    auto& edge_dbound) {
                    auto& edge_internal = edge_dbound.edge_data.edge_internal;

                    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

                    uint ndof_local  = edge_dbound.boundary.data.get_ndof();
                    uint ndof_global = edge_dbound.edge_data.get_ndof();

                    uint local_offset  = boundary.local_dof_offset;
                    uint global_offset = edge_internal.global_dof_offset;

                    subvector(rhs_global_temp, global_offset * SWE::n_variables, ndof_global * SWE::n_variables) =
                        edge_internal.rhs_global;

                    for (uint i = 0; i < ndof_local * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_local.add_triplet(local_offset * SWE::n_variables + i,
                                                               global_offset * SWE::n_variables + j,
                                                               boundary.delta_hat_local(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_local * SWE::n_variables; ++j) {
                            sparse_delta_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            local_offset * SWE::n_variables + j,
                                                            boundary.delta_global(i, j));
                        }
                    }
                });

            sparse_delta_local_inv.get_sparse_matrix(delta_local_inv_temp);
            sparse_delta_hat_local.get_sparse_matrix(delta_hat_local_temp);
            sparse_delta_global.get_sparse_matrix(delta_global_temp);
            sparse_delta_hat_global.get_sparse_matrix(delta_hat_global_temp);

            delta_local_inv += delta_local_inv_temp;
            delta_hat_local += delta_hat_local_temp;
            rhs_local += rhs_local_temp;

            delta_global += delta_global_temp;
            delta_hat_global += delta_hat_global_temp;
            rhs_global += rhs_global_temp;
        }

        delta_hat_global -= delta_global * delta_local_inv * delta_hat_local;

        rhs_global -= delta_global * delta_local_inv * rhs_local;

        solve_sle(delta_hat_global, rhs_global);

        rhs_local = delta_local_inv * (rhs_local - delta_hat_local * rhs_global);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([&rhs_local, &sim_units, su_id](auto& elt) {
                const uint stage = sim_units[su_id]->stepper.GetStage();

                auto& state = elt.data.state[stage + 1];

                uint ndof         = elt.data.get_ndof();
                uint local_offset = elt.data.internal.local_dof_offset;

                auto sol_local = subvector(rhs_local, local_offset * SWE::n_variables, ndof * SWE::n_variables);

                for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
                    state.q(SWE::Variables::ze, dof) += sol_local[3 * dof];
                    state.q(SWE::Variables::qx, dof) += sol_local[3 * dof + 1];
                    state.q(SWE::Variables::qy, dof) += sol_local[3 * dof + 2];
                }
            });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global](auto& edge_int) {
                auto& edge_state = edge_int.edge_data.edge_state;

                uint ndof          = edge_int.edge_data.get_ndof();
                uint global_offset = edge_int.edge_data.edge_internal.global_dof_offset;

                auto sol_global = subvector(rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

                for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
                    edge_state.q_hat(SWE::Variables::ze, dof) += sol_global[3 * dof];
                    edge_state.q_hat(SWE::Variables::qx, dof) += sol_global[3 * dof + 1];
                    edge_state.q_hat(SWE::Variables::qy, dof) += sol_global[3 * dof + 2];
                }
            });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global](auto& edge_bound) {
                auto& edge_state = edge_bound.edge_data.edge_state;

                uint ndof          = edge_bound.edge_data.get_ndof();
                uint global_offset = edge_bound.edge_data.edge_internal.global_dof_offset;

                auto sol_global = subvector(rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

                for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
                    edge_state.q_hat(SWE::Variables::ze, dof) += sol_global[3 * dof];
                    edge_state.q_hat(SWE::Variables::qx, dof) += sol_global[3 * dof + 1];
                    edge_state.q_hat(SWE::Variables::qy, dof) += sol_global[3 * dof + 2];
                }
            });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&rhs_global](auto& edge_dbound) {
                auto& edge_state = edge_dbound.edge_data.edge_state;

                uint ndof          = edge_dbound.edge_data.get_ndof();
                uint global_offset = edge_dbound.edge_data.edge_internal.global_dof_offset;

                auto sol_global = subvector(rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

                for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
                    edge_state.q_hat(SWE::Variables::ze, dof) += sol_global[3 * dof];
                    edge_state.q_hat(SWE::Variables::qx, dof) += sol_global[3 * dof + 1];
                    edge_state.q_hat(SWE::Variables::qy, dof) += sol_global[3 * dof + 2];
                }
            });
        }
    }

#pragma omp barrier

    return false;
}
}
}

#endif