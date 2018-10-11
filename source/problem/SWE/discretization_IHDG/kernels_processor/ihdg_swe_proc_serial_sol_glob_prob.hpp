#ifndef IHDG_SWE_PROC_SERIAL_SOL_GLOB_PROB_HPP
#define IHDG_SWE_PROC_SERIAL_SOL_GLOB_PROB_HPP

namespace SWE {
namespace IHDG {
bool Problem::serial_solve_global_problem(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    auto& global_data = discretization.global_data;

    SparseMatrix<double>& delta_hat_global = global_data.delta_hat_global;
    DynVector<double>& rhs_global          = global_data.rhs_global;

    SparseMatrixMeta<double> sparse_delta_hat_global;

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        internal.delta_local_inv = inverse(internal.delta_local);
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global, &sparse_delta_hat_global](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        edge_internal.delta_hat_global -=
            boundary_in.delta_global * internal_in.delta_local_inv * boundary_in.delta_hat_local +
            boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_ex.delta_hat_local;

        edge_internal.rhs_global -= boundary_in.delta_global * internal_in.delta_local_inv * internal_in.rhs_local +
                                    boundary_ex.delta_global * internal_ex.delta_local_inv * internal_ex.rhs_local;

        subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size()) = edge_internal.rhs_global;

        for (uint i = 0; i < global_dof_indx.size(); ++i) {
            for (uint j = 0; j < global_dof_indx.size(); ++j) {
                sparse_delta_hat_global.add_triplet(
                    global_dof_indx[i], global_dof_indx[j], edge_internal.delta_hat_global(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_in)
                continue;

            auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary_in.delta_global * internal_in.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }

        for (uint bound_id = 0; bound_id < edge_int.interface.data_ex.get_nbound(); ++bound_id) {
            if (bound_id == edge_int.interface.bound_id_ex)
                continue;

            auto& boundary_con = edge_int.interface.data_ex.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary_ex.delta_global * internal_ex.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global, &sparse_delta_hat_global](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        edge_internal.delta_hat_global -= boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

        edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

        subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size()) = edge_internal.rhs_global;

        for (uint i = 0; i < global_dof_indx.size(); ++i) {
            for (uint j = 0; j < global_dof_indx.size(); ++j) {
                sparse_delta_hat_global.add_triplet(
                    global_dof_indx[i], global_dof_indx[j], edge_internal.delta_hat_global(i, j));
            }
        }

        for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
            if (bound_id == edge_bound.boundary.bound_id)
                continue;

            auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

            edge_internal.delta_hat_global =
                -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

            std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_dof_indx[i], global_dof_con_indx[j], edge_internal.delta_hat_global(i, j));
                }
            }
        }
    });

    sparse_delta_hat_global.get_sparse_matrix(delta_hat_global);

    solve_sle(delta_hat_global, rhs_global);

    discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        auto del_q_hat = subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size());

        internal_in.rhs_local -= boundary_in.delta_hat_local * del_q_hat;
        internal_ex.rhs_local -= boundary_ex.delta_hat_local * del_q_hat;

        edge_state.q_hat +=
            reshape<double, SWE::n_variables, SO::ColumnMajor>(del_q_hat, edge_int.edge_data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

        auto del_q_hat = subvector(rhs_global, (uint)global_dof_indx[0], (uint)global_dof_indx.size());

        internal.rhs_local -= boundary.delta_hat_local * del_q_hat;

        edge_state.q_hat +=
            reshape<double, SWE::n_variables, SO::ColumnMajor>(del_q_hat, edge_bound.edge_data.get_ndof());
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        auto& internal = elt.data.internal;

        internal.rhs_local = internal.delta_local_inv * internal.rhs_local;

        state.q += reshape<double, SWE::n_variables, SO::ColumnMajor>(internal.rhs_local, elt.data.get_ndof());
    });

    double delta_norm = norm(rhs_global) / rows(delta_hat_global);

    if (delta_norm < 1e-8) {
        return true;
    }

    return false;
}
}
}

#endif