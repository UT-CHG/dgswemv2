#ifndef IHDG_SWE_PROC_SOLVE_GLOBAL_PROB_HPP
#define IHDG_SWE_PROC_SOLVE_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename SimulationType>
bool Problem::solve_global_problem(SimulationType* simulation) {
    SparseMatrixData<double> sparse_delta_local_inv;
    SparseMatrixData<double> sparse_delta_hat_local;
    SparseMatrixData<double> sparse_delta_global;
    SparseMatrixData<double> sparse_delta_hat_global;

    simulation->mesh.CallForEachElement([simulation, &sparse_delta_local_inv](auto& elt) {
        auto& internal = elt.data.internal;

        uint ndof         = elt.data.get_ndof();
        uint local_offset = internal.local_dof_offset;

        subvector(simulation->rhs_local, local_offset, ndof * SWE::n_variables) = internal.rhs_local;

        internal.delta_local_inv = inverse(internal.delta_local);

        for (uint i = 0; i < ndof * SWE::n_variables; ++i) {
            for (uint j = 0; j < ndof * SWE::n_variables; ++j) {
                sparse_delta_local_inv.add_triplet(local_offset + i, local_offset + j, internal.delta_local_inv(i, j));
            }
        }
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface(
        [simulation, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](auto& edge_int) {
            auto& edge_internal = edge_int.edge_data.edge_internal;

            auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
            auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

            uint ndof_local_in = edge_int.interface.data_in.get_ndof();
            uint ndof_local_ex = edge_int.interface.data_ex.get_ndof();
            uint ndof_global   = edge_int.edge_data.get_ndof();

            uint local_offset_in = boundary_in.local_dof_offset;
            uint local_offset_ex = boundary_ex.local_dof_offset;
            uint global_offset   = edge_internal.global_dof_offset;

            subvector(simulation->rhs_global, global_offset, ndof_global * SWE::n_variables) = edge_internal.rhs_global;

            for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_offset + i, global_offset + j, edge_internal.delta_hat_global(i, j));
                }
            }

            for (uint i = 0; i < ndof_local_in * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                    sparse_delta_hat_local.add_triplet(
                        local_offset_in + i, global_offset + j, boundary_in.delta_hat_local(i, j));
                }
            }

            for (uint i = 0; i < ndof_local_ex * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                    sparse_delta_hat_local.add_triplet(
                        local_offset_ex + i, global_offset + j, boundary_ex.delta_hat_local(i, j));
                }
            }

            for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_local_in * SWE::n_variables; ++j) {
                    sparse_delta_global.add_triplet(
                        global_offset + i, local_offset_in + j, boundary_in.delta_global(i, j));
                }
            }

            for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_local_ex * SWE::n_variables; ++j) {
                    sparse_delta_global.add_triplet(
                        global_offset + i, local_offset_ex + j, boundary_ex.delta_global(i, j));
                }
            }
        });

    simulation->mesh_skeleton.CallForEachEdgeBoundary(
        [simulation, &sparse_delta_hat_local, &sparse_delta_global, &sparse_delta_hat_global](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            uint ndof_local  = edge_bound.boundary.data.get_ndof();
            uint ndof_global = edge_bound.edge_data.get_ndof();

            uint local_offset  = boundary.local_dof_offset;
            uint global_offset = edge_internal.global_dof_offset;

            subvector(simulation->rhs_global, global_offset, ndof_global * SWE::n_variables) = edge_internal.rhs_global;

            for (uint i = 0; i < ndof_local * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                    sparse_delta_hat_local.add_triplet(
                        local_offset + i, global_offset + j, boundary.delta_hat_local(i, j));
                }
            }

            for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_global * SWE::n_variables; ++j) {
                    sparse_delta_hat_global.add_triplet(
                        global_offset + i, global_offset + j, edge_internal.delta_hat_global(i, j));
                }
            }

            for (uint i = 0; i < ndof_global * SWE::n_variables; ++i) {
                for (uint j = 0; j < ndof_local * SWE::n_variables; ++j) {
                    sparse_delta_global.add_triplet(global_offset + i, local_offset + j, boundary.delta_global(i, j));
                }
            }
        });

    simulation->delta_local_inv.setFromTriplets(sparse_delta_local_inv.data.begin(), sparse_delta_local_inv.data.end());
    simulation->delta_hat_local.setFromTriplets(sparse_delta_hat_local.data.begin(), sparse_delta_hat_local.data.end());

    simulation->delta_global.setFromTriplets(sparse_delta_global.data.begin(), sparse_delta_global.data.end());
    simulation->delta_hat_global.setFromTriplets(sparse_delta_hat_global.data.begin(),
                                                 sparse_delta_hat_global.data.end());

    simulation->delta_hat_global = simulation->delta_hat_global -
                                   simulation->delta_global * simulation->delta_local_inv * simulation->delta_hat_local;

    simulation->rhs_global =
        simulation->rhs_global - simulation->delta_global * simulation->delta_local_inv * simulation->rhs_local;

    solve_sle(simulation->delta_hat_global, simulation->rhs_global);

    simulation->rhs_local =
        simulation->delta_local_inv * (simulation->rhs_local - simulation->delta_hat_local * simulation->rhs_global);

    simulation->mesh.CallForEachElement([simulation](auto& elt) {
        const uint stage = simulation->stepper.GetStage();

        auto& state = elt.data.state[stage + 1];

        uint ndof         = elt.data.get_ndof();
        uint local_offset = elt.data.internal.local_dof_offset;

        auto rhs_local = subvector(simulation->rhs_local, local_offset, ndof * SWE::n_variables);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q(SWE::Variables::ze, dof) += rhs_local[3 * dof];
            state.q(SWE::Variables::qx, dof) += rhs_local[3 * dof + 1];
            state.q(SWE::Variables::qy, dof) += rhs_local[3 * dof + 2];
        }
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface([simulation](auto& edge_int) {
        auto& edge_state = edge_int.edge_data.edge_state;

        uint ndof = edge_int.edge_data.get_ndof();
        uint global_offset = edge_int.edge_data.edge_internal.global_dof_offset;

        auto rhs_global = subvector(simulation->rhs_global, global_offset, ndof * SWE::n_variables);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) += rhs_global[3 * dof];
            edge_state.q_hat(SWE::Variables::qx, dof) += rhs_global[3 * dof + 1];
            edge_state.q_hat(SWE::Variables::qy, dof) += rhs_global[3 * dof + 2];
        }
    });

    simulation->mesh_skeleton.CallForEachEdgeBoundary([simulation](auto& edge_bound) {
        auto& edge_state = edge_bound.edge_data.edge_state;

        uint ndof = edge_bound.edge_data.get_ndof();
        uint global_offset = edge_bound.edge_data.edge_internal.global_dof_offset;

        auto rhs_global = subvector(simulation->rhs_global, global_offset, ndof * SWE::n_variables);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat(SWE::Variables::ze, dof) += rhs_global[3 * dof];
            edge_state.q_hat(SWE::Variables::qx, dof) += rhs_global[3 * dof + 1];
            edge_state.q_hat(SWE::Variables::qy, dof) += rhs_global[3 * dof + 2];
        }
    });

    double delta_norm = norm(simulation->rhs_global) / simulation->rhs_global.size();

    if (delta_norm < 1e-8) {
        return true;
    }

    return false;
}
}
}

#endif