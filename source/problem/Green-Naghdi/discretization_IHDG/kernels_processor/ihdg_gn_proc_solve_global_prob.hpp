#ifndef IHDG_GN_PROC_SOLVE_GLOBAL_PROB_HPP
#define IHDG_GN_PROC_SOLVE_GLOBAL_PROB_HPP

namespace GN {
namespace IHDG {
void Problem::solve_global_dc_problem(const RKStepper& stepper, HDGDiscretization<Problem>& discretization) {
    SparseMatrixMeta<double> sparse_w1_w1_hat;

    SparseMatrixMeta<double> sparse_w2_w1;
    SparseMatrixMeta<double> sparse_w2_w2_inv;
    SparseMatrixMeta<double> sparse_w2_w1_hat;

    SparseMatrixMeta<double> sparse_w1_hat_w1;
    SparseMatrixMeta<double> sparse_w1_hat_w2;
    SparseMatrixMeta<double> sparse_w1_hat_w1_hat;

    discretization.mesh.CallForEachElement([&discretization, &sparse_w2_w1, &sparse_w2_w2_inv](auto& elt) {
        auto& internal = elt.data.internal;

        uint ndof         = elt.data.get_ndof();
        uint local_offset = internal.local_dof_offset;

        internal.w2_w2_inv = inverse(internal.w2_w2);

        internal.w1_w1 -= internal.w1_w2 * internal.w2_w2_inv * internal.w2_w1;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].w1_w1_hat -=
                internal.w1_w2 * internal.w2_w2_inv * elt.data.boundary[bound_id].w2_w1_hat;

            solve_sle(internal.w1_w1, elt.data.boundary[bound_id].w1_w1_hat);
        }

        solve_sle(internal.w1_w1, internal.w1_rhs);

        subvector(discretization.global_data.w1_rhs, local_offset * GN::n_dimensions, ndof * GN::n_dimensions) =
            internal.w1_rhs;

        for (uint i = 0; i < ndof; ++i) {
            for (uint j = 0; j < ndof * GN::n_dimensions; ++j) {
                sparse_w2_w1.add_triplet(local_offset + i, local_offset * GN::n_dimensions + j, internal.w2_w1(i, j));
            }
        }

        for (uint i = 0; i < ndof; ++i) {
            for (uint j = 0; j < ndof; ++j) {
                sparse_w2_w2_inv.add_triplet(local_offset + i, local_offset + j, internal.w2_w2_inv(i, j));
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&discretization,
                                                           &sparse_w1_w1_hat,
                                                           &sparse_w2_w1_hat,
                                                           &sparse_w1_hat_w1,
                                                           &sparse_w1_hat_w2,
                                                           &sparse_w1_hat_w1_hat](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint ndof_local_in = edge_int.interface.data_in.get_ndof();
        uint ndof_local_ex = edge_int.interface.data_ex.get_ndof();
        uint ndof_global   = edge_int.edge_data.get_ndof();

        uint local_offset_in = boundary_in.local_dof_offset;
        uint local_offset_ex = boundary_ex.local_dof_offset;
        uint global_offset   = edge_internal.global_dof_offset;

        for (uint i = 0; i < ndof_global * GN::n_dimensions; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w1_hat_w1_hat.add_triplet(global_offset * GN::n_dimensions + i,
                                                 global_offset * GN::n_dimensions + j,
                                                 edge_internal.w1_hat_w1_hat(i, j));
            }
        }

        for (uint i = 0; i < ndof_local_in * GN::n_dimensions; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w1_w1_hat.add_triplet(local_offset_in * GN::n_dimensions + i,
                                             global_offset * GN::n_dimensions + j,
                                             boundary_in.w1_w1_hat(i, j));

                sparse_w1_hat_w1.add_triplet(global_offset * GN::n_dimensions + j,
                                             local_offset_in * GN::n_dimensions + i,
                                             boundary_in.w1_hat_w1(j, i));
            }
        }

        for (uint i = 0; i < ndof_local_in; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w2_w1_hat.add_triplet(
                    local_offset_in + i, global_offset * GN::n_dimensions + j, boundary_in.w2_w1_hat(i, j));

                sparse_w1_hat_w2.add_triplet(
                    global_offset * GN::n_dimensions + j, local_offset_in + i, boundary_in.w1_hat_w2(j, i));
            }
        }

        for (uint i = 0; i < ndof_local_ex * GN::n_dimensions; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w1_w1_hat.add_triplet(local_offset_ex * GN::n_dimensions + i,
                                             global_offset * GN::n_dimensions + j,
                                             boundary_ex.w1_w1_hat(i, j));

                sparse_w1_hat_w1.add_triplet(global_offset * GN::n_dimensions + j,
                                             local_offset_ex * GN::n_dimensions + i,
                                             boundary_ex.w1_hat_w1(j, i));
            }
        }

        for (uint i = 0; i < ndof_local_ex; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w2_w1_hat.add_triplet(
                    local_offset_ex + i, global_offset * GN::n_dimensions + j, boundary_ex.w2_w1_hat(i, j));

                sparse_w1_hat_w2.add_triplet(
                    global_offset * GN::n_dimensions + j, local_offset_ex + i, boundary_ex.w1_hat_w2(j, i));
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&discretization,
                                                          &sparse_w1_w1_hat,
                                                          &sparse_w2_w1_hat,
                                                          &sparse_w1_hat_w1,
                                                          &sparse_w1_hat_w2,
                                                          &sparse_w1_hat_w1_hat](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint ndof_local  = edge_bound.boundary.data.get_ndof();
        uint ndof_global = edge_bound.edge_data.get_ndof();

        uint local_offset  = boundary.local_dof_offset;
        uint global_offset = edge_internal.global_dof_offset;

        for (uint i = 0; i < ndof_global * GN::n_dimensions; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w1_hat_w1_hat.add_triplet(global_offset * GN::n_dimensions + i,
                                                 global_offset * GN::n_dimensions + j,
                                                 edge_internal.w1_hat_w1_hat(i, j));
            }
        }

        for (uint i = 0; i < ndof_local * GN::n_dimensions; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w1_w1_hat.add_triplet(local_offset * GN::n_dimensions + i,
                                             global_offset * GN::n_dimensions + j,
                                             boundary.w1_w1_hat(i, j));

                sparse_w1_hat_w1.add_triplet(global_offset * GN::n_dimensions + j,
                                             local_offset * GN::n_dimensions + i,
                                             boundary.w1_hat_w1(j, i));
            }
        }

        for (uint i = 0; i < ndof_local; ++i) {
            for (uint j = 0; j < ndof_global * GN::n_dimensions; ++j) {
                sparse_w2_w1_hat.add_triplet(
                    local_offset + i, global_offset * GN::n_dimensions + j, boundary.w2_w1_hat(i, j));

                // sparse_w1_hat_w2.add_triplet(
                //    global_offset * GN::n_dimensions + j, local_offset + i, boundary.w1_hat_w2(j, i));
            }
        }
    });

    sparse_w1_w1_hat.get_sparse_matrix(discretization.global_data.w1_w1_hat);

    sparse_w2_w1.get_sparse_matrix(discretization.global_data.w2_w1);
    sparse_w2_w2_inv.get_sparse_matrix(discretization.global_data.w2_w2_inv);
    sparse_w2_w1_hat.get_sparse_matrix(discretization.global_data.w2_w1_hat);

    sparse_w1_hat_w1.get_sparse_matrix(discretization.global_data.w1_hat_w1);
    sparse_w1_hat_w2.get_sparse_matrix(discretization.global_data.w1_hat_w2);
    sparse_w1_hat_w1_hat.get_sparse_matrix(discretization.global_data.w1_hat_w1_hat);

    discretization.global_data.w1_hat_w1 -=
        discretization.global_data.w1_hat_w2 * discretization.global_data.w2_w2_inv * discretization.global_data.w2_w1;
    discretization.global_data.w1_hat_w1_hat -= discretization.global_data.w1_hat_w2 *
                                                discretization.global_data.w2_w2_inv *
                                                discretization.global_data.w2_w1_hat;

    discretization.global_data.w1_hat_w1_hat -=
        discretization.global_data.w1_hat_w1 * discretization.global_data.w1_w1_hat;
    discretization.global_data.w1_hat_rhs = -discretization.global_data.w1_hat_w1 * discretization.global_data.w1_rhs;

    solve_sle(discretization.global_data.w1_hat_w1_hat, discretization.global_data.w1_hat_rhs);

    discretization.global_data.w1_rhs -= discretization.global_data.w1_w1_hat * discretization.global_data.w1_hat_rhs;

    discretization.mesh.CallForEachElement([&discretization, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        uint ndof         = elt.data.get_ndof();
        uint local_offset = elt.data.internal.local_dof_offset;

        auto w1_rhs =
            subvector(discretization.global_data.w1_rhs, local_offset * GN::n_dimensions, ndof * GN::n_dimensions);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.w1(GlobalCoord::x, dof) = w1_rhs[2 * dof];
            state.w1(GlobalCoord::y, dof) = w1_rhs[2 * dof + 1];
        }
    });
}
}
}

#endif