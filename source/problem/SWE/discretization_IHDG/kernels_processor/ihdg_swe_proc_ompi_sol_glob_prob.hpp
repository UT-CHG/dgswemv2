#ifndef IHDG_SWE_PROC_OMPI_SOL_GLOB_PROB_HPP
#define IHDG_SWE_PROC_OMPI_SOL_GLOB_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename OMPISimUnitType>
bool Problem::ompi_solve_global_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    uint global_dof_offset = 2 * 3136;

    SparseMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;

    delta_hat_global.resize(global_dof_offset * SWE::n_variables, global_dof_offset * SWE::n_variables);
    rhs_global.resize(global_dof_offset * SWE::n_variables);

    set_constant(rhs_global, 0.0);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        SparseMatrix<double> delta_hat_global_temp;
        DynVector<double> rhs_global_temp;

        delta_hat_global_temp.resize(global_dof_offset * SWE::n_variables, global_dof_offset * SWE::n_variables);
        rhs_global_temp.resize(global_dof_offset * SWE::n_variables);

        set_constant(rhs_global_temp, 0.0);

        SparseMatrixMeta<double> sparse_delta_hat_global;

        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& internal = elt.data.internal;

            internal.delta_local_inv = inverse(internal.delta_local);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
            [&rhs_global_temp, &sparse_delta_hat_global](auto& edge_int) {
                auto& edge_internal = edge_int.edge_data.edge_internal;

                auto& internal_in = edge_int.interface.data_in.internal;
                auto& internal_ex = edge_int.interface.data_ex.internal;

                auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
                auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

                uint ndof_global   = edge_int.edge_data.get_ndof();
                uint global_offset = edge_internal.global_dof_offset;

                edge_internal.delta_hat_global -=
                    boundary_in.delta_global * internal_in.delta_local_inv * boundary_in.delta_hat_local +
                    boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_ex.delta_hat_local;

                edge_internal.rhs_global -=
                    boundary_in.delta_global * internal_in.delta_local_inv * internal_in.rhs_local +
                    boundary_ex.delta_global * internal_ex.delta_local_inv * internal_ex.rhs_local;

                subvector(rhs_global_temp, global_offset * SWE::n_variables, ndof_global * SWE::n_variables) =
                    edge_internal.rhs_global;

                for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                    for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                        sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            global_offset * SWE::n_variables + j,
                                                            edge_internal.delta_hat_global(i, j));
                    }
                }

                for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
                    if (bound_id == edge_int.interface.bound_id_in)
                        continue;

                    auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

                    edge_internal.delta_hat_global =
                        -boundary_in.delta_global * internal_in.delta_local_inv * boundary_con.delta_hat_local;

                    uint global_offset_con = boundary_con.global_dof_offset;

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset_con * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }
                }

                for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
                    if (bound_id == edge_int.interface.bound_id_ex)
                        continue;

                    auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];

                    edge_internal.delta_hat_global =
                        -boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_con.delta_hat_local;

                    uint global_offset_con = boundary_con.global_dof_offset;

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset_con * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }
                }
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
            [&rhs_global_temp, &sparse_delta_hat_global](auto& edge_bound) {
                auto& edge_internal = edge_bound.edge_data.edge_internal;

                auto& internal = edge_bound.boundary.data.internal;
                auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

                uint ndof_global   = edge_bound.edge_data.get_ndof();
                uint global_offset = edge_internal.global_dof_offset;

                edge_internal.delta_hat_global -=
                    boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

                edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

                subvector(rhs_global_temp, global_offset * SWE::n_variables, ndof_global * SWE::n_variables) =
                    edge_internal.rhs_global;

                for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                    for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                        sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            global_offset * SWE::n_variables + j,
                                                            edge_internal.delta_hat_global(i, j));
                    }
                }

                for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
                    if (bound_id == edge_bound.boundary.bound_id)
                        continue;

                    auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

                    edge_internal.delta_hat_global =
                        -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

                    uint global_offset_con = boundary_con.global_dof_offset;

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset_con * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }
                }
            });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
            [&rhs_global_temp, &sparse_delta_hat_global](auto& edge_dbound) {
                auto& edge_internal = edge_dbound.edge_data.edge_internal;

                auto& internal = edge_dbound.boundary.data.internal;
                auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

                uint ndof_global   = edge_dbound.edge_data.get_ndof();
                uint global_offset = edge_internal.global_dof_offset;

                edge_internal.delta_hat_global -=
                    boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

                edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

                subvector(rhs_global_temp, global_offset * SWE::n_variables, ndof_global * SWE::n_variables) =
                    edge_internal.rhs_global;

                for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                    for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                        sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                            global_offset * SWE::n_variables + j,
                                                            edge_internal.delta_hat_global(i, j));
                    }
                }

                for (uint bound_id = 0; bound_id < edge_dbound.boundary.data.get_nbound(); ++bound_id) {
                    if (bound_id == edge_dbound.boundary.bound_id)
                        continue;

                    auto& boundary_con = edge_dbound.boundary.data.boundary[bound_id];

                    edge_internal.delta_hat_global =
                        -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

                    uint global_offset_con = boundary_con.global_dof_offset;

                    for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                        for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                            sparse_delta_hat_global.add_triplet(global_offset * SWE::n_variables + i,
                                                                global_offset_con * SWE::n_variables + j,
                                                                edge_internal.delta_hat_global(i, j));
                        }
                    }
                }
            });

        sparse_delta_hat_global.get_sparse_matrix(delta_hat_global_temp);

        delta_hat_global += delta_hat_global_temp;
        rhs_global += rhs_global_temp;
    }

    solve_sle(delta_hat_global, rhs_global);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global](auto& edge_int) {
            auto& edge_state = edge_int.edge_data.edge_state;

            auto& internal_in = edge_int.interface.data_in.internal;
            auto& internal_ex = edge_int.interface.data_ex.internal;

            auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
            auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

            uint ndof          = edge_int.edge_data.get_ndof();
            uint global_offset = edge_int.edge_data.edge_internal.global_dof_offset;

            auto del_q_hat = subvector(rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

            internal_in.rhs_local -= boundary_in.delta_hat_local * del_q_hat;
            internal_ex.rhs_local -= boundary_ex.delta_hat_local * del_q_hat;

            for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
                edge_state.q_hat(SWE::Variables::ze, dof) += del_q_hat[3 * dof];
                edge_state.q_hat(SWE::Variables::qx, dof) += del_q_hat[3 * dof + 1];
                edge_state.q_hat(SWE::Variables::qy, dof) += del_q_hat[3 * dof + 2];
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global](auto& edge_bound) {
            auto& edge_state = edge_bound.edge_data.edge_state;

            auto& internal = edge_bound.boundary.data.internal;
            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            uint ndof          = edge_bound.edge_data.get_ndof();
            uint global_offset = edge_bound.edge_data.edge_internal.global_dof_offset;

            auto del_q_hat = subvector(rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

            internal.rhs_local -= boundary.delta_hat_local * del_q_hat;

            for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
                edge_state.q_hat(SWE::Variables::ze, dof) += del_q_hat[3 * dof];
                edge_state.q_hat(SWE::Variables::qx, dof) += del_q_hat[3 * dof + 1];
                edge_state.q_hat(SWE::Variables::qy, dof) += del_q_hat[3 * dof + 2];
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&rhs_global](auto& edge_dbound) {
            auto& edge_state = edge_dbound.edge_data.edge_state;

            auto& internal = edge_dbound.boundary.data.internal;
            auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            uint ndof          = edge_dbound.edge_data.get_ndof();
            uint global_offset = edge_dbound.edge_data.edge_internal.global_dof_offset;

            auto del_q_hat = subvector(rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

            internal.rhs_local -= boundary.delta_hat_local * del_q_hat;

            for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
                edge_state.q_hat(SWE::Variables::ze, dof) += del_q_hat[3 * dof];
                edge_state.q_hat(SWE::Variables::qx, dof) += del_q_hat[3 * dof + 1];
                edge_state.q_hat(SWE::Variables::qy, dof) += del_q_hat[3 * dof + 2];
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state = elt.data.state[stage + 1];

            auto& internal = elt.data.internal;

            internal.rhs_local *= internal.delta_local_inv;

            for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
                state.q(SWE::Variables::ze, dof) += internal.rhs_local[3 * dof];
                state.q(SWE::Variables::qx, dof) += internal.rhs_local[3 * dof + 1];
                state.q(SWE::Variables::qy, dof) += internal.rhs_local[3 * dof + 2];
            }
        });
    }

    double delta_norm = norm(rhs_global) / rhs_global.size();

    if (delta_norm < 1e-8) {
        return true;
    }

    return false;
}
}
}

#endif