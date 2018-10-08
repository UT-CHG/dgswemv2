#ifndef EHDG_GN_PROC_SERIAL_SOL_GLOB_PROB_HPP
#define EHDG_GN_PROC_SERIAL_SOL_GLOB_PROB_HPP

namespace GN {
namespace EHDG {
void Problem::solve_global_dc_problem(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    auto& global_data = discretization.global_data;

    SparseMatrix<double>& w1_hat_w1_hat = global_data.w1_hat_w1_hat;
    DynVector<double>& w1_hat_rhs       = global_data.w1_hat_rhs;

    SparseMatrixMeta<double> sparse_w1_hat_w1_hat;

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

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

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        std::vector<uint>& dc_global_dof_indx = edge_internal.dc_global_dof_indx;

        boundary_in.w1_hat_w1 -= boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * internal_in.w2_w1;

        boundary_ex.w1_hat_w1 -= boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * internal_ex.w2_w1;

        edge_internal.w1_hat_w1_hat -= boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_in.w2_w1_hat +
                                       boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * boundary_ex.w2_w1_hat +
                                       boundary_in.w1_hat_w1 * boundary_in.w1_w1_hat +
                                       boundary_ex.w1_hat_w1 * boundary_ex.w1_w1_hat;

        subvector(w1_hat_rhs, (uint)dc_global_dof_indx[0], (uint)dc_global_dof_indx.size()) =
            -(boundary_in.w1_hat_w1 * internal_in.w1_rhs + boundary_ex.w1_hat_w1 * internal_ex.w1_rhs);

        for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
            for (uint j = 0; j < dc_global_dof_indx.size(); ++j) {
                sparse_w1_hat_w1_hat.add_triplet(
                    dc_global_dof_indx[i], dc_global_dof_indx[j], edge_internal.w1_hat_w1_hat(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_in)
                continue;

            auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

            edge_internal.w1_hat_w1_hat = -(boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_con.w2_w1_hat +
                                            boundary_in.w1_hat_w1 * boundary_con.w1_w1_hat);

            std::vector<uint>& dc_global_dof_con_indx = boundary_con.dc_global_dof_indx;

            for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                for (uint j = 0; j < dc_global_dof_con_indx.size(); ++j) {
                    sparse_w1_hat_w1_hat.add_triplet(
                        dc_global_dof_indx[i], dc_global_dof_con_indx[j], edge_internal.w1_hat_w1_hat(i, j));
                }
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_ex)
                continue;

            auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];

            edge_internal.w1_hat_w1_hat = -(boundary_ex.w1_hat_w2 * internal_ex.w2_w2_inv * boundary_con.w2_w1_hat +
                                            boundary_ex.w1_hat_w1 * boundary_con.w1_w1_hat);

            std::vector<uint>& dc_global_dof_con_indx = boundary_con.dc_global_dof_indx;

            for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                for (uint j = 0; j < dc_global_dof_con_indx.size(); ++j) {
                    sparse_w1_hat_w1_hat.add_triplet(
                        dc_global_dof_indx[i], dc_global_dof_con_indx[j], edge_internal.w1_hat_w1_hat(i, j));
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&w1_hat_rhs, &sparse_w1_hat_w1_hat](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        std::vector<uint>& dc_global_dof_indx = edge_internal.dc_global_dof_indx;

        /* boundary.w1_hat_w1 -= boundary.w1_hat_w2 * internal.w2_w2_inv * internal.w2_w1; */

        edge_internal.w1_hat_w1_hat -=
            /* boundary.w1_hat_w2 * internal.w2_w2_inv * boundary.w2_w1_hat + */ boundary.w1_hat_w1 *
            boundary.w1_w1_hat;

        subvector(w1_hat_rhs, (uint)dc_global_dof_indx[0], (uint)dc_global_dof_indx.size()) =
            -boundary.w1_hat_w1 * internal.w1_rhs;

        for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
            for (uint j = 0; j < dc_global_dof_indx.size(); ++j) {
                sparse_w1_hat_w1_hat.add_triplet(
                    dc_global_dof_indx[i], dc_global_dof_indx[j], edge_internal.w1_hat_w1_hat(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
            if (bound_id == edge_bound.boundary.bound_id)
                continue;

            auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

            edge_internal.w1_hat_w1_hat = -(/* boundary.w1_hat_w2 * internal.w2_w2_inv * boundary_con.w2_w1_hat + */
                                            boundary.w1_hat_w1 * boundary_con.w1_w1_hat);

            std::vector<uint>& dc_global_dof_con_indx = boundary_con.dc_global_dof_indx;

            for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                for (uint j = 0; j < dc_global_dof_con_indx.size(); ++j) {
                    sparse_w1_hat_w1_hat.add_triplet(
                        dc_global_dof_indx[i], dc_global_dof_con_indx[j], edge_internal.w1_hat_w1_hat(i, j));
                }
            }
        }
    });

    sparse_w1_hat_w1_hat.get_sparse_matrix(w1_hat_w1_hat);

    solve_sle(w1_hat_w1_hat, w1_hat_rhs);

    discretization.mesh.CallForEachElement([&w1_hat_rhs, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            std::vector<uint>& dc_global_dof_indx = elt.data.boundary[bound_id].dc_global_dof_indx;

            auto w1_hat = subvector(w1_hat_rhs, (uint)dc_global_dof_indx[0], (uint)dc_global_dof_indx.size());

            internal.w1_rhs -= elt.data.boundary[bound_id].w1_w1_hat * w1_hat;
        }

        state.w1 = reshape<double, GN::n_dimensions, SO::ColumnMajor>(internal.w1_rhs, elt.data.get_ndof());
    });
}
}
}

#endif