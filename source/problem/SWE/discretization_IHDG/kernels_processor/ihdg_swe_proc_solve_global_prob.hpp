#ifndef IHDG_SWE_PROC_SOLVE_GLOBAL_PROB_HPP
#define IHDG_SWE_PROC_SOLVE_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
bool Problem::solve_global_problem(const RKStepper& stepper, HDGDiscretization<Problem>& discretization) {
    SparseMatrixMeta<double> sparse_delta_local_inv;
    SparseMatrixMeta<double> sparse_delta_hat_local;
    SparseMatrixMeta<double> sparse_delta_global;
    SparseMatrixMeta<double> sparse_delta_hat_global;

    discretization.mesh.CallForEachElement([&discretization, &sparse_delta_local_inv](auto& elt) {
        auto& internal = elt.data.internal;

        uint ndof         = elt.data.get_ndof();
        uint local_offset = internal.local_dof_offset;

        subvector(discretization.global_data.rhs_local, local_offset * SWE::n_variables, ndof * SWE::n_variables) =
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

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&discretization, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](auto& edge_int) {
            auto& edge_internal = edge_int.edge_data.edge_internal;

            auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
            auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

            uint ndof_local_in = edge_int.interface.data_in.get_ndof();
            uint ndof_local_ex = edge_int.interface.data_ex.get_ndof();
            uint ndof_global   = edge_int.edge_data.get_ndof();

            uint local_offset_in = boundary_in.local_dof_offset;
            uint local_offset_ex = boundary_ex.local_dof_offset;
            uint global_offset   = edge_internal.global_dof_offset;

            subvector(discretization.global_data.rhs_global,
                      global_offset * SWE::n_variables,
                      ndof_global * SWE::n_variables) = edge_internal.rhs_global;

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

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&discretization, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            uint ndof_local  = edge_bound.boundary.data.get_ndof();
            uint ndof_global = edge_bound.edge_data.get_ndof();

            uint local_offset  = boundary.local_dof_offset;
            uint global_offset = edge_internal.global_dof_offset;

            subvector(discretization.global_data.rhs_global,
                      global_offset * SWE::n_variables,
                      ndof_global * SWE::n_variables) = edge_internal.rhs_global;

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

    sparse_delta_local_inv.get_sparse_matrix(discretization.global_data.delta_local_inv);
    sparse_delta_hat_local.get_sparse_matrix(discretization.global_data.delta_hat_local);
    sparse_delta_global.get_sparse_matrix(discretization.global_data.delta_global);
    sparse_delta_hat_global.get_sparse_matrix(discretization.global_data.delta_hat_global);

    discretization.global_data.delta_hat_global -= discretization.global_data.delta_global *
                                                   discretization.global_data.delta_local_inv *
                                                   discretization.global_data.delta_hat_local;

    discretization.global_data.rhs_global -= discretization.global_data.delta_global *
                                             discretization.global_data.delta_local_inv *
                                             discretization.global_data.rhs_local;

    solve_sle(discretization.global_data.delta_hat_global, discretization.global_data.rhs_global);

    discretization.global_data.rhs_local =
        discretization.global_data.delta_local_inv *
        (discretization.global_data.rhs_local -
         discretization.global_data.delta_hat_local * discretization.global_data.rhs_global);

    discretization.mesh.CallForEachElement([&discretization, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        uint ndof         = elt.data.get_ndof();
        uint local_offset = elt.data.internal.local_dof_offset;

        auto rhs_local =
            subvector(discretization.global_data.rhs_local, local_offset * SWE::n_variables, ndof * SWE::n_variables);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q(SWE::Variables::ze, dof) += rhs_local[3 * dof];
            state.q(SWE::Variables::qx, dof) += rhs_local[3 * dof + 1];
            state.q(SWE::Variables::qy, dof) += rhs_local[3 * dof + 2];
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&discretization](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        uint ndof          = edge_int.edge_data.get_ndof();
        uint global_offset = edge_int.edge_data.edge_internal.global_dof_offset;

        auto rhs_global =
            subvector(discretization.global_data.rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) += rhs_global[3 * dof];
            edge_state.q_hat(SWE::Variables::qx, dof) += rhs_global[3 * dof + 1];
            edge_state.q_hat(SWE::Variables::qy, dof) += rhs_global[3 * dof + 2];
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&discretization](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        uint ndof          = edge_bound.edge_data.get_ndof();
        uint global_offset = edge_bound.edge_data.edge_internal.global_dof_offset;

        auto rhs_global =
            subvector(discretization.global_data.rhs_global, global_offset * SWE::n_variables, ndof * SWE::n_variables);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) += rhs_global[3 * dof];
            edge_state.q_hat(SWE::Variables::qx, dof) += rhs_global[3 * dof + 1];
            edge_state.q_hat(SWE::Variables::qy, dof) += rhs_global[3 * dof + 2];
        }
    });

    double delta_norm = norm(discretization.global_data.rhs_global) / discretization.global_data.rhs_global.size();

    if (delta_norm < 1e-8) {
        return true;
    }

    return false;
}
}
}

#endif