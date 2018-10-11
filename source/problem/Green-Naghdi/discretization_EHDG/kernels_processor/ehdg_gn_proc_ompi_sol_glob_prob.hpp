#ifndef EHDG_GN_PROC_OMPI_SOL_GLOB_PROB_HPP
#define EHDG_GN_PROC_OMPI_SOL_GLOB_PROB_HPP

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::ompi_solve_global_dc_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    auto& global_data = sim_units[0]->discretization.global_data;

    Mat& w1_hat_w1_hat = global_data.w1_hat_w1_hat;
    Vec& w1_hat_rhs    = global_data.w1_hat_rhs;

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
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
    }

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&w1_hat_rhs,
                                                                                 &w1_hat_w1_hat](auto& edge_int) {
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

            edge_internal.w1_hat_rhs =
                -(boundary_in.w1_hat_w1 * internal_in.w1_rhs + boundary_ex.w1_hat_w1 * internal_ex.w1_rhs);

            std::vector<double> w1_hat_rhs_container(dc_global_dof_indx.size());
            std::vector<double> w1_hat_w1_hat_container(dc_global_dof_indx.size() * dc_global_dof_indx.size());

            for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                w1_hat_rhs_container[i] = edge_internal.w1_hat_rhs[i];

                for (uint j = 0; j < dc_global_dof_indx.size(); ++j) {
                    w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                }
            }

            VecSetValues(w1_hat_rhs,
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         &w1_hat_rhs_container.front(),
                         ADD_VALUES);

            MatSetValues(w1_hat_w1_hat,
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         &w1_hat_w1_hat_container.front(),
                         ADD_VALUES);

            for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
                if (bound_id == edge_int.interface.bound_id_in)
                    continue;

                auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

                edge_internal.w1_hat_w1_hat = -(boundary_in.w1_hat_w2 * internal_in.w2_w2_inv * boundary_con.w2_w1_hat +
                                                boundary_in.w1_hat_w1 * boundary_con.w1_w1_hat);

                std::vector<uint>& dc_global_dof_con_indx = boundary_con.dc_global_dof_indx;

                for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                    for (uint j = 0; j < dc_global_dof_con_indx.size(); ++j) {
                        w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                    }
                }

                MatSetValues(w1_hat_w1_hat,
                             dc_global_dof_indx.size(),
                             (int*)&dc_global_dof_indx.front(),
                             dc_global_dof_con_indx.size(),
                             (int*)&dc_global_dof_con_indx.front(),
                             &w1_hat_w1_hat_container.front(),
                             ADD_VALUES);
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
                        w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                    }
                }

                MatSetValues(w1_hat_w1_hat,
                             dc_global_dof_indx.size(),
                             (int*)&dc_global_dof_indx.front(),
                             dc_global_dof_con_indx.size(),
                             (int*)&dc_global_dof_con_indx.front(),
                             &w1_hat_w1_hat_container.front(),
                             ADD_VALUES);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&w1_hat_rhs,
                                                                                &w1_hat_w1_hat](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            auto& internal = edge_bound.boundary.data.internal;
            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            std::vector<uint>& dc_global_dof_indx = edge_internal.dc_global_dof_indx;

            /* boundary.w1_hat_w1 -= boundary.w1_hat_w2 * internal.w2_w2_inv * internal.w2_w1; */

            edge_internal.w1_hat_w1_hat -=
                /* boundary.w1_hat_w2 * internal.w2_w2_inv * boundary.w2_w1_hat + */ boundary.w1_hat_w1 *
                boundary.w1_w1_hat;

            edge_internal.w1_hat_rhs = -boundary.w1_hat_w1 * internal.w1_rhs;

            std::vector<double> w1_hat_rhs_container(dc_global_dof_indx.size());
            std::vector<double> w1_hat_w1_hat_container(dc_global_dof_indx.size() * dc_global_dof_indx.size());

            for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                w1_hat_rhs_container[i] = edge_internal.w1_hat_rhs[i];

                for (uint j = 0; j < dc_global_dof_indx.size(); ++j) {
                    w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                }
            }

            VecSetValues(w1_hat_rhs,
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         &w1_hat_rhs_container.front(),
                         ADD_VALUES);

            MatSetValues(w1_hat_w1_hat,
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         &w1_hat_w1_hat_container.front(),
                         ADD_VALUES);

            for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
                if (bound_id == edge_bound.boundary.bound_id)
                    continue;

                auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

                edge_internal.w1_hat_w1_hat = -(/* boundary.w1_hat_w2 * internal.w2_w2_inv * boundary_con.w2_w1_hat + */
                                                boundary.w1_hat_w1 * boundary_con.w1_w1_hat);

                std::vector<uint>& dc_global_dof_con_indx = boundary_con.dc_global_dof_indx;

                for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                    for (uint j = 0; j < dc_global_dof_con_indx.size(); ++j) {
                        w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                    }
                }

                MatSetValues(w1_hat_w1_hat,
                             dc_global_dof_indx.size(),
                             (int*)&dc_global_dof_indx.front(),
                             dc_global_dof_con_indx.size(),
                             (int*)&dc_global_dof_con_indx.front(),
                             &w1_hat_w1_hat_container.front(),
                             ADD_VALUES);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&w1_hat_rhs,
                                                                                   &w1_hat_w1_hat](auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;

            auto& internal = edge_dbound.boundary.data.internal;
            auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            std::vector<uint>& dc_global_dof_indx = edge_internal.dc_global_dof_indx;

            boundary.w1_hat_w1 -= boundary.w1_hat_w2 * internal.w2_w2_inv * internal.w2_w1;

            edge_internal.w1_hat_w1_hat -=
                boundary.w1_hat_w2 * internal.w2_w2_inv * boundary.w2_w1_hat + boundary.w1_hat_w1 * boundary.w1_w1_hat;

            edge_internal.w1_hat_rhs = -boundary.w1_hat_w1 * internal.w1_rhs;

            std::vector<double> w1_hat_rhs_container(dc_global_dof_indx.size());
            std::vector<double> w1_hat_w1_hat_container(dc_global_dof_indx.size() * dc_global_dof_indx.size());

            for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                w1_hat_rhs_container[i] = edge_internal.w1_hat_rhs[i];

                for (uint j = 0; j < dc_global_dof_indx.size(); ++j) {
                    w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                }
            }

            VecSetValues(w1_hat_rhs,
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         &w1_hat_rhs_container.front(),
                         ADD_VALUES);

            MatSetValues(w1_hat_w1_hat,
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         dc_global_dof_indx.size(),
                         (int*)&dc_global_dof_indx.front(),
                         &w1_hat_w1_hat_container.front(),
                         ADD_VALUES);

            for (uint bound_id = 0; bound_id < edge_dbound.boundary.data.get_nbound(); ++bound_id) {
                if (bound_id == edge_dbound.boundary.bound_id)
                    continue;

                auto& boundary_con = edge_dbound.boundary.data.boundary[bound_id];

                edge_internal.w1_hat_w1_hat = -(boundary.w1_hat_w2 * internal.w2_w2_inv * boundary_con.w2_w1_hat +
                                                boundary.w1_hat_w1 * boundary_con.w1_w1_hat);

                std::vector<uint>& dc_global_dof_con_indx = boundary_con.dc_global_dof_indx;

                for (uint i = 0; i < dc_global_dof_indx.size(); ++i) {
                    for (uint j = 0; j < dc_global_dof_con_indx.size(); ++j) {
                        w1_hat_w1_hat_container[i * dc_global_dof_indx.size() + j] = edge_internal.w1_hat_w1_hat(i, j);
                    }
                }

                MatSetValues(w1_hat_w1_hat,
                             dc_global_dof_indx.size(),
                             (int*)&dc_global_dof_indx.front(),
                             dc_global_dof_con_indx.size(),
                             (int*)&dc_global_dof_con_indx.front(),
                             &w1_hat_w1_hat_container.front(),
                             ADD_VALUES);
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

    KSPSolve(dc_ksp, w1_hat_rhs, w1_hat_rhs);

    VecScatterBegin(dc_scatter, w1_hat_rhs, dc_sol, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(dc_scatter, w1_hat_rhs, dc_sol, INSERT_VALUES, SCATTER_FORWARD);

    double* sol_ptr;
    VecGetArray(dc_sol, &sol_ptr);

    int sol_size;
    VecGetLocalSize(dc_sol, &sol_size);

    auto solution = vector_from_array(sol_ptr, sol_size);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&solution](auto& edge_int) {
            auto& edge_internal = edge_int.edge_data.edge_internal;

            auto& internal_in = edge_int.interface.data_in.internal;
            auto& internal_ex = edge_int.interface.data_ex.internal;

            auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
            auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

            uint n_global_dofs = edge_internal.dc_global_dof_indx.size();

            auto w1_hat = subvector(solution, edge_internal.dc_sol_offset, n_global_dofs);

            internal_in.w1_rhs -= boundary_in.w1_w1_hat * w1_hat;
            internal_ex.w1_rhs -= boundary_ex.w1_w1_hat * w1_hat;
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&solution](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            auto& internal = edge_bound.boundary.data.internal;
            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            uint n_global_dofs = edge_internal.dc_global_dof_indx.size();

            auto w1_hat = subvector(solution, edge_internal.dc_sol_offset, n_global_dofs);

            internal.w1_rhs -= boundary.w1_w1_hat * w1_hat;
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&solution](auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;

            auto& internal = edge_dbound.boundary.data.internal;
            auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            uint n_global_dofs = edge_internal.dc_global_dof_indx.size();

            auto w1_hat = subvector(solution, edge_internal.dc_sol_offset, n_global_dofs);

            internal.w1_rhs -= boundary.w1_w1_hat * w1_hat;
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = elt.data.state[stage];
            auto& internal = elt.data.internal;

            state.w1 = reshape<double, GN::n_dimensions, SO::ColumnMajor>(internal.w1_rhs, elt.data.get_ndof());
        });
    }

    VecRestoreArray(dc_sol, &sol_ptr);

    MatZeroEntries(w1_hat_w1_hat);
    VecZeroEntries(w1_hat_rhs);
}
}
}

#endif