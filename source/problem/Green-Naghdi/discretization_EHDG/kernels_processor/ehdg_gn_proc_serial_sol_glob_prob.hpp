#ifndef EHDG_GN_PROC_SERIAL_SOL_GLOB_PROB_HPP
#define EHDG_GN_PROC_SERIAL_SOL_GLOB_PROB_HPP

namespace GN {
namespace EHDG {
void Problem::serial_solve_global_dc_problem(ProblemDiscretizationType& discretization,
                                             ProblemGlobalDataType& global_data,
                                             const ESSPRKStepper& stepper) {
    SparseMatrix<double>& w1_hat_w1_hat = global_data.w1_hat_w1_hat;
    DynVector<double>& w1_hat_rhs       = global_data.w1_hat_rhs;
    SparseMatrixMeta<double> sparse_w1_hat_w1_hat;

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal     = elt.data.internal;
        internal.w2_w2_inv = inverse(internal.w2_w2);
        internal.w1_w1 -= internal.w1_w2 * internal.w2_w2_inv * internal.w2_w1;
        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].w1_w1_hat -=
                internal.w1_w2 * internal.w2_w2_inv * elt.data.boundary[bound_id].w2_w1_hat;
            solve_sle(internal.w1_w1, elt.data.boundary[bound_id].w1_w1_hat);
        }
        solve_sle(internal.w1_w1, internal.w1_rhs);
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&w1_hat_rhs, &sparse_w1_hat_w1_hat](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;
        auto& internal_in   = edge_int.interface.data_in.internal;
        auto& internal_ex   = edge_int.interface.data_ex.internal;
        auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        boundary_in.w1_hat_w1 -= boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * internal_in.w2_w1;
        boundary_ex.w1_hat_w1 -= boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * internal_ex.w2_w1;
        edge_internal.w1_hat_w1_hat -= boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_in.w2_w1_hat +
                                       boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * boundary_ex.w2_w1_hat +
                                       boundary_in.w1_hat_w1 * boundary_in.w1_w1_hat +
                                       boundary_ex.w1_hat_w1 * boundary_ex.w1_w1_hat;

        const uint n_global_dofs = edge_int.edge_data.get_ndof() * GN::n_dimensions;
        subvector(w1_hat_rhs, edge_internal.dc_global_dof_indx * n_global_dofs, n_global_dofs) =
            -(boundary_in.w1_hat_w1 * internal_in.w1_rhs + boundary_ex.w1_hat_w1 * internal_ex.w1_rhs);
        for (uint i = 0; i < n_global_dofs; ++i) {
            for (uint j = 0; j < n_global_dofs; ++j) {
                sparse_w1_hat_w1_hat.add_triplet(edge_internal.dc_global_dof_indx * n_global_dofs + i,
                                                 edge_internal.dc_global_dof_indx * n_global_dofs + j,
                                                 edge_internal.w1_hat_w1_hat(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_in)
                continue;
            auto& boundary_con          = edge_int.interface.data_in.boundary[bound_id];
            edge_internal.w1_hat_w1_hat = -(boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_con.w2_w1_hat +
                                            boundary_in.w1_hat_w1 * boundary_con.w1_w1_hat);
            for (uint i = 0; i < n_global_dofs; ++i) {
                for (uint j = 0; j < n_global_dofs; ++j) {
                    sparse_w1_hat_w1_hat.add_triplet(edge_internal.dc_global_dof_indx * n_global_dofs + i,
                                                     boundary_con.dc_global_dof_indx * n_global_dofs + j,
                                                     edge_internal.w1_hat_w1_hat(i, j));
                }
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_ex)
                continue;
            auto& boundary_con          = edge_int.interface.data_ex.boundary[bound_id];
            edge_internal.w1_hat_w1_hat = -(boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * boundary_con.w2_w1_hat +
                                            boundary_ex.w1_hat_w1 * boundary_con.w1_w1_hat);
            for (uint i = 0; i < n_global_dofs; ++i) {
                for (uint j = 0; j < n_global_dofs; ++j) {
                    sparse_w1_hat_w1_hat.add_triplet(edge_internal.dc_global_dof_indx * n_global_dofs + i,
                                                     boundary_con.dc_global_dof_indx * n_global_dofs + j,
                                                     edge_internal.w1_hat_w1_hat(i, j));
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&w1_hat_rhs, &sparse_w1_hat_w1_hat](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;
        auto& internal      = edge_bound.boundary.data.internal;
        auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        /* boundary.w1_hat_w1 -= boundary.w1_hat_w2 * internal.w2_w2_inv * internal.w2_w1; */
        edge_internal.w1_hat_w1_hat -=
            /* boundary.w1_hat_w2 * internal.w2_w2_inv * boundary.w2_w1_hat + */ boundary.w1_hat_w1 *
            boundary.w1_w1_hat;

        const uint n_global_dofs = edge_bound.edge_data.get_ndof() * GN::n_dimensions;
        subvector(w1_hat_rhs, edge_internal.dc_global_dof_indx * n_global_dofs, n_global_dofs) =
            -boundary.w1_hat_w1 * internal.w1_rhs;
        for (uint i = 0; i < n_global_dofs; ++i) {
            for (uint j = 0; j < n_global_dofs; ++j) {
                sparse_w1_hat_w1_hat.add_triplet(edge_internal.dc_global_dof_indx * n_global_dofs + i,
                                                 edge_internal.dc_global_dof_indx * n_global_dofs + j,
                                                 edge_internal.w1_hat_w1_hat(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
            if (bound_id == edge_bound.boundary.bound_id)
                continue;
            auto& boundary_con          = edge_bound.boundary.data.boundary[bound_id];
            edge_internal.w1_hat_w1_hat = -(/* boundary.w1_hat_w2 * internal.w2_w2_inv * boundary_con.w2_w1_hat + */
                                            boundary.w1_hat_w1 * boundary_con.w1_w1_hat);
            for (uint i = 0; i < n_global_dofs; ++i) {
                for (uint j = 0; j < n_global_dofs; ++j) {
                    sparse_w1_hat_w1_hat.add_triplet(edge_internal.dc_global_dof_indx * n_global_dofs + i,
                                                     boundary_con.dc_global_dof_indx * n_global_dofs + j,
                                                     edge_internal.w1_hat_w1_hat(i, j));
                }
            }
        }
    });

    sparse_w1_hat_w1_hat.get_sparse_matrix(w1_hat_w1_hat);

    solve_sle(w1_hat_w1_hat, w1_hat_rhs);

    discretization.mesh_skeleton.CallForEachEdgeInterface([&w1_hat_rhs](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;
        auto& internal_in   = edge_int.interface.data_in.internal;
        auto& internal_ex   = edge_int.interface.data_ex.internal;
        auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        const uint n_global_dofs = edge_int.edge_data.get_ndof() * GN::n_dimensions;
        const auto w1_hat = subvector(w1_hat_rhs, edge_internal.dc_local_dof_indx * n_global_dofs, n_global_dofs);
        internal_in.w1_rhs -= boundary_in.w1_w1_hat * w1_hat;
        internal_ex.w1_rhs -= boundary_ex.w1_w1_hat * w1_hat;
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&w1_hat_rhs](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;
        auto& internal      = edge_bound.boundary.data.internal;
        auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        const uint n_global_dofs = edge_bound.edge_data.get_ndof() * GN::n_dimensions;
        const auto w1_hat = subvector(w1_hat_rhs, edge_internal.dc_local_dof_indx * n_global_dofs, n_global_dofs);
        internal.w1_rhs -= boundary.w1_w1_hat * w1_hat;
    });

    discretization.mesh.CallForEachElement([&w1_hat_rhs, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& internal   = elt.data.internal;
        state.w1         = reshape<double, GN::n_dimensions, SO::ColumnMajor>(internal.w1_rhs, elt.data.get_ndof());
    });
}
}
}

#endif