#ifndef IHDG_SWE_PROC_OMPI_SOL_GLOB_PROB_HPP
#define IHDG_SWE_PROC_OMPI_SOL_GLOB_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename OMPISimUnitType>
bool Problem::ompi_solve_global_problem(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
    uint n_global_dofs = sim_units[0]->discretization.n_global_dofs;

    Mat delta_hat_global;
    MatCreate(MPI_COMM_SELF, &delta_hat_global);
    MatSetSizes(delta_hat_global, PETSC_DECIDE, PETSC_DECIDE, n_global_dofs, n_global_dofs);
    MatSetUp(delta_hat_global);

    Vec rhs_global;
    VecCreateMPI(MPI_COMM_SELF, PETSC_DECIDE, n_global_dofs, &rhs_global);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& internal = elt.data.internal;

            internal.delta_local_inv = inverse(internal.delta_local);
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&rhs_global,
                                                                                 &delta_hat_global](auto& edge_int) {
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

            std::vector<double> rhs_global_container(global_dof_indx.size());
            std::vector<double> delta_hat_global_container(global_dof_indx.size() * global_dof_indx.size());

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                rhs_global_container[i] = edge_internal.rhs_global[i];

                for (uint j = 0; j < global_dof_indx.size(); ++j) {
                    delta_hat_global_container[i * global_dof_indx.size() + j] = edge_internal.delta_hat_global(i, j);
                }
            }

            VecSetValues(rhs_global,
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         &rhs_global_container.front(),
                         ADD_VALUES);

            MatSetValues(delta_hat_global,
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         &delta_hat_global_container.front(),
                         ADD_VALUES);

            for (uint bound_id = 0; bound_id < edge_int.interface.data_in.get_nbound(); ++bound_id) {
                if (bound_id == edge_int.interface.bound_id_in)
                    continue;

                auto& boundary_con = edge_int.interface.data_in.boundary[bound_id];

                edge_internal.delta_hat_global =
                    -boundary_in.delta_global * internal_in.delta_local_inv * boundary_con.delta_hat_local;

                std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

                for (uint i = 0; i < global_dof_indx.size(); ++i) {
                    for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                        delta_hat_global_container[i * global_dof_indx.size() + j] =
                            edge_internal.delta_hat_global(i, j);
                    }
                }

                MatSetValues(delta_hat_global,
                             global_dof_indx.size(),
                             (int*)&global_dof_indx.front(),
                             global_dof_con_indx.size(),
                             (int*)&global_dof_con_indx.front(),
                             &delta_hat_global_container.front(),
                             ADD_VALUES);
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
                        delta_hat_global_container[i * global_dof_indx.size() + j] =
                            edge_internal.delta_hat_global(i, j);
                    }
                }

                MatSetValues(delta_hat_global,
                             global_dof_indx.size(),
                             (int*)&global_dof_indx.front(),
                             global_dof_con_indx.size(),
                             (int*)&global_dof_con_indx.front(),
                             &delta_hat_global_container.front(),
                             ADD_VALUES);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&rhs_global,
                                                                                &delta_hat_global](auto& edge_bound) {
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            auto& internal = edge_bound.boundary.data.internal;
            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

            edge_internal.delta_hat_global -=
                boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

            edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

            std::vector<double> rhs_global_container(global_dof_indx.size());
            std::vector<double> delta_hat_global_container(global_dof_indx.size() * global_dof_indx.size());

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                rhs_global_container[i] = edge_internal.rhs_global[i];

                for (uint j = 0; j < global_dof_indx.size(); ++j) {
                    delta_hat_global_container[i * global_dof_indx.size() + j] = edge_internal.delta_hat_global(i, j);
                }
            }

            VecSetValues(rhs_global,
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         &rhs_global_container.front(),
                         ADD_VALUES);

            MatSetValues(delta_hat_global,
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         &delta_hat_global_container.front(),
                         ADD_VALUES);

            for (uint bound_id = 0; bound_id < edge_bound.boundary.data.get_nbound(); ++bound_id) {
                if (bound_id == edge_bound.boundary.bound_id)
                    continue;

                auto& boundary_con = edge_bound.boundary.data.boundary[bound_id];

                edge_internal.delta_hat_global =
                    -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

                std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

                for (uint i = 0; i < global_dof_indx.size(); ++i) {
                    for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                        delta_hat_global_container[i * global_dof_indx.size() + j] =
                            edge_internal.delta_hat_global(i, j);
                    }
                }

                MatSetValues(delta_hat_global,
                             global_dof_indx.size(),
                             (int*)&global_dof_indx.front(),
                             global_dof_con_indx.size(),
                             (int*)&global_dof_con_indx.front(),
                             &delta_hat_global_container.front(),
                             ADD_VALUES);
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&rhs_global, &delta_hat_global](
                                                                                      auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;

            auto& internal = edge_dbound.boundary.data.internal;
            auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

            edge_internal.delta_hat_global -=
                boundary.delta_global * internal.delta_local_inv * boundary.delta_hat_local;

            edge_internal.rhs_global -= boundary.delta_global * internal.delta_local_inv * internal.rhs_local;

            std::vector<double> rhs_global_container(global_dof_indx.size());
            std::vector<double> delta_hat_global_container(global_dof_indx.size() * global_dof_indx.size());

            for (uint i = 0; i < global_dof_indx.size(); ++i) {
                rhs_global_container[i] = edge_internal.rhs_global[i];

                for (uint j = 0; j < global_dof_indx.size(); ++j) {
                    delta_hat_global_container[i * global_dof_indx.size() + j] = edge_internal.delta_hat_global(i, j);
                }
            }

            VecSetValues(rhs_global,
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         &rhs_global_container.front(),
                         ADD_VALUES);

            MatSetValues(delta_hat_global,
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         global_dof_indx.size(),
                         (int*)&global_dof_indx.front(),
                         &delta_hat_global_container.front(),
                         ADD_VALUES);

            for (uint bound_id = 0; bound_id < edge_dbound.boundary.data.get_nbound(); ++bound_id) {
                if (bound_id == edge_dbound.boundary.bound_id)
                    continue;

                auto& boundary_con = edge_dbound.boundary.data.boundary[bound_id];

                edge_internal.delta_hat_global =
                    -boundary.delta_global * internal.delta_local_inv * boundary_con.delta_hat_local;

                std::vector<uint>& global_dof_con_indx = boundary_con.global_dof_indx;

                for (uint i = 0; i < global_dof_indx.size(); ++i) {
                    for (uint j = 0; j < global_dof_con_indx.size(); ++j) {
                        delta_hat_global_container[i * global_dof_indx.size() + j] =
                            edge_internal.delta_hat_global(i, j);
                    }
                }

                MatSetValues(delta_hat_global,
                             global_dof_indx.size(),
                             (int*)&global_dof_indx.front(),
                             global_dof_con_indx.size(),
                             (int*)&global_dof_con_indx.front(),
                             &delta_hat_global_container.front(),
                             ADD_VALUES);
            }
        });
    }

    MatAssemblyBegin(delta_hat_global, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(delta_hat_global, MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(rhs_global);
    VecAssemblyEnd(rhs_global);

    KSP ksp;
    PC pc;

    KSPCreate(MPI_COMM_SELF, &ksp);
    KSPSetOperators(ksp, delta_hat_global, delta_hat_global);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);

    KSPSolve(ksp, rhs_global, rhs_global);

    double* sol;
    VecGetArray(rhs_global, &sol);

    Eigen::Map<DynVector<double>> solution(sol, n_global_dofs);

    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface([&solution](auto& edge_int) {
            auto& edge_state    = edge_int.edge_data.edge_state;
            auto& edge_internal = edge_int.edge_data.edge_internal;

            auto& internal_in = edge_int.interface.data_in.internal;
            auto& internal_ex = edge_int.interface.data_ex.internal;

            auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
            auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

            std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

            auto del_q_hat = subvector(solution, global_dof_indx[0], global_dof_indx.size());

            internal_in.rhs_local -= boundary_in.delta_hat_local * del_q_hat;
            internal_ex.rhs_local -= boundary_ex.delta_hat_local * del_q_hat;

            for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
                edge_state.q_hat(SWE::Variables::ze, dof) += del_q_hat[3 * dof];
                edge_state.q_hat(SWE::Variables::qx, dof) += del_q_hat[3 * dof + 1];
                edge_state.q_hat(SWE::Variables::qy, dof) += del_q_hat[3 * dof + 2];
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary([&solution](auto& edge_bound) {
            auto& edge_state    = edge_bound.edge_data.edge_state;
            auto& edge_internal = edge_bound.edge_data.edge_internal;

            auto& internal = edge_bound.boundary.data.internal;
            auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

            std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

            auto del_q_hat = subvector(solution, global_dof_indx[0], global_dof_indx.size());

            internal.rhs_local -= boundary.delta_hat_local * del_q_hat;

            for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
                edge_state.q_hat(SWE::Variables::ze, dof) += del_q_hat[3 * dof];
                edge_state.q_hat(SWE::Variables::qx, dof) += del_q_hat[3 * dof + 1];
                edge_state.q_hat(SWE::Variables::qy, dof) += del_q_hat[3 * dof + 2];
            }
        });

        sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed([&solution](auto& edge_dbound) {
            auto& edge_state    = edge_dbound.edge_data.edge_state;
            auto& edge_internal = edge_dbound.edge_data.edge_internal;

            auto& internal = edge_dbound.boundary.data.internal;
            auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            std::vector<uint>& global_dof_indx = edge_internal.global_dof_indx;

            auto del_q_hat = subvector(solution, global_dof_indx[0], global_dof_indx.size());

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

            internal.rhs_local = internal.delta_local_inv * internal.rhs_local;

            for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
                state.q(SWE::Variables::ze, dof) += internal.rhs_local[3 * dof];
                state.q(SWE::Variables::qx, dof) += internal.rhs_local[3 * dof + 1];
                state.q(SWE::Variables::qy, dof) += internal.rhs_local[3 * dof + 2];
            }
        });
    }

    double delta_norm = norm(solution) / n_global_dofs;

    VecRestoreArray(rhs_global, &sol);

    MatDestroy(&delta_hat_global);
    VecDestroy(&rhs_global);
    KSPDestroy(&ksp);

    if (delta_norm < 1e-8) {
        return true;
    }

    return false;
}
}
}

#endif