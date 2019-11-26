#ifndef RECON_AVG_DBATH_OMPI_HPP
#define RECON_AVG_DBATH_OMPI_HPP

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_dbath(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, ProblemGlobalDataType& global_data) {
#pragma omp barrier
#pragma omp master
    {
        std::set<uint> nodeIDs;
        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&nodeIDs](auto& elt) { nodeIDs.insert(elt.GetNodeID().begin(), elt.GetNodeID().end()); });
        }
        global_data.local_nodeIDs = std::vector<int>(nodeIDs.begin(), nodeIDs.end());
        std::map<uint, uint> nodeIDmap;
        uint it = 0;
        for (auto nodeID : nodeIDs) {
            nodeIDmap[nodeID] = it++;
        }

        uint local_max_nodeID  = *std::max_element(nodeIDs.begin(), nodeIDs.end());
        uint global_max_nodeID = 0;
        MPI_Allreduce(&local_max_nodeID, &global_max_nodeID, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

        VecCreateMPI(MPI_COMM_WORLD,
                     (global_max_nodeID + 1) * GN::n_dddbath_terms,
                     PETSC_DECIDE,
                     &(global_data.global_derivatives_at_node));
        VecSetBlockSize(global_data.global_derivatives_at_node, GN::n_dddbath_terms);
        VecCreateSeq(MPI_COMM_SELF,
                     global_data.local_nodeIDs.size() * GN::n_dddbath_terms,
                     &(global_data.local_derivatives_at_node));
        global_data.derivatives_at_node = DynVector<double>(global_data.local_nodeIDs.size() * GN::n_dddbath_terms);

        ISCreateBlock(MPI_COMM_SELF,
                      GN::n_dddbath_terms,
                      global_data.local_nodeIDs.size(),
                      (int*)&global_data.local_nodeIDs.front(),
                      PETSC_COPY_VALUES,
                      &(global_data.derivatives_from));
        ISCreateStride(
            MPI_COMM_SELF, global_data.local_nodeIDs.size() * GN::n_dddbath_terms, 0, 1, &(global_data.derivatives_to));
        VecScatterCreate(global_data.global_derivatives_at_node,
                         global_data.derivatives_from,
                         global_data.local_derivatives_at_node,
                         global_data.derivatives_to,
                         &(global_data.derivatives_scatter));

        VecCreateMPI(MPI_COMM_WORLD, global_max_nodeID + 1, PETSC_DECIDE, &(global_data.global_node_mult));
        VecCreateSeq(MPI_COMM_SELF, global_data.local_nodeIDs.size(), &(global_data.local_node_mult));

        ISCreateGeneral(MPI_COMM_SELF,
                        global_data.local_nodeIDs.size(),
                        (int*)&global_data.local_nodeIDs.front(),
                        PETSC_COPY_VALUES,
                        &(global_data.node_mult_from));
        ISCreateStride(MPI_COMM_SELF, global_data.local_nodeIDs.size(), 0, 1, &(global_data.node_mult_to));
        VecScatterCreate(global_data.global_node_mult,
                         global_data.node_mult_from,
                         global_data.local_node_mult,
                         global_data.node_mult_to,
                         &(global_data.node_mult_scatter));

        VecSet(global_data.global_derivatives_at_node, 0.0);
        VecSet(global_data.global_node_mult, 0.0);
        set_constant(global_data.derivatives_at_node, 0.0);

        DynVector<double> node_mult(global_data.local_nodeIDs.size());
        set_constant(node_mult, 0.0);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([&global_data, &node_mult, &nodeIDmap](auto& elt) {
                auto& derivative = elt.data.derivative;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    derivative.local_nodeID[node] = nodeIDmap[elt.GetNodeID()[node]];
                    subvector(global_data.derivatives_at_node,
                              derivative.local_nodeID[node] * GN::n_dddbath_terms,
                              GN::n_dimensions) += derivative.dbath_at_baryctr;
                    node_mult[derivative.local_nodeID[node]] += 1.0;
                }
            });
        }

        VecSetValuesBlocked(global_data.global_derivatives_at_node,
                            global_data.local_nodeIDs.size(),
                            &global_data.local_nodeIDs.front(),
                            global_data.derivatives_at_node.data(),
                            ADD_VALUES);
        VecSetValues(global_data.global_node_mult,
                     global_data.local_nodeIDs.size(),
                     &global_data.local_nodeIDs.front(),
                     node_mult.data(),
                     ADD_VALUES);

        VecAssemblyBegin(global_data.global_derivatives_at_node);
        VecAssemblyBegin(global_data.global_node_mult);
        VecAssemblyEnd(global_data.global_derivatives_at_node);
        VecAssemblyEnd(global_data.global_node_mult);

        VecScatterBegin(global_data.derivatives_scatter,
                        global_data.global_derivatives_at_node,
                        global_data.local_derivatives_at_node,
                        INSERT_VALUES,
                        SCATTER_FORWARD);
        VecScatterBegin(global_data.node_mult_scatter,
                        global_data.global_node_mult,
                        global_data.local_node_mult,
                        INSERT_VALUES,
                        SCATTER_FORWARD);
        VecScatterEnd(global_data.derivatives_scatter,
                      global_data.global_derivatives_at_node,
                      global_data.local_derivatives_at_node,
                      INSERT_VALUES,
                      SCATTER_FORWARD);
        VecScatterEnd(global_data.node_mult_scatter,
                      global_data.global_node_mult,
                      global_data.local_node_mult,
                      INSERT_VALUES,
                      SCATTER_FORWARD);

        double* d_ptr;
        VecGetArray(global_data.local_derivatives_at_node, &d_ptr);

        double* m_ptr;
        VecGetArray(global_data.local_node_mult, &m_ptr);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([d_ptr, m_ptr](auto& elt) {
                auto& state      = elt.data.state[0];
                auto& derivative = elt.data.derivative;
                auto& internal   = elt.data.internal;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    derivative.node_mult[node] = m_ptr[derivative.local_nodeID[node]];
                    for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                        derivative.dbath_lin(dbath, node) =
                            d_ptr[GN::n_dddbath_terms * derivative.local_nodeID[node] + dbath] /
                            derivative.node_mult[node];
                    }
                }
                state.dbath          = elt.ProjectLinearToBasis(derivative.dbath_lin);
                internal.dbath_at_gp = elt.ComputeUgp(state.dbath);
            });
        }

        VecRestoreArray(global_data.local_derivatives_at_node, &d_ptr);
        VecRestoreArray(global_data.local_node_mult, &m_ptr);
    }
#pragma omp barrier
}

template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_ddbath(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, ProblemGlobalDataType& global_data) {
#pragma omp barrier
#pragma omp master
    {
        VecSet(global_data.global_derivatives_at_node, 0.0);
        set_constant(global_data.derivatives_at_node, 0.0);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([&global_data](auto& elt) {
                auto& derivative = elt.data.derivative;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    subvector(global_data.derivatives_at_node,
                              derivative.local_nodeID[node] * GN::n_dddbath_terms,
                              GN::n_ddbath_terms) += derivative.ddbath_at_baryctr;
                }
            });
        }

        VecSetValuesBlocked(global_data.global_derivatives_at_node,
                            global_data.local_nodeIDs.size(),
                            &global_data.local_nodeIDs.front(),
                            global_data.derivatives_at_node.data(),
                            ADD_VALUES);

        VecAssemblyBegin(global_data.global_derivatives_at_node);
        VecAssemblyEnd(global_data.global_derivatives_at_node);

        VecScatterBegin(global_data.derivatives_scatter,
                        global_data.global_derivatives_at_node,
                        global_data.local_derivatives_at_node,
                        INSERT_VALUES,
                        SCATTER_FORWARD);
        VecScatterEnd(global_data.derivatives_scatter,
                      global_data.global_derivatives_at_node,
                      global_data.local_derivatives_at_node,
                      INSERT_VALUES,
                      SCATTER_FORWARD);

        double* d_ptr;
        VecGetArray(global_data.local_derivatives_at_node, &d_ptr);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([d_ptr](auto& elt) {
                auto& state      = elt.data.state[0];
                auto& derivative = elt.data.derivative;
                auto& internal   = elt.data.internal;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                        derivative.ddbath_lin(ddbath, node) =
                            d_ptr[GN::n_dddbath_terms * derivative.local_nodeID[node] + ddbath] /
                            derivative.node_mult[node];
                    }
                }
                state.ddbath          = elt.ProjectLinearToBasis(derivative.ddbath_lin);
                internal.ddbath_at_gp = elt.ComputeUgp(state.ddbath);
            });
        }

        VecRestoreArray(global_data.local_derivatives_at_node, &d_ptr);
    }
#pragma omp barrier
}

template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_dddbath(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, ProblemGlobalDataType& global_data) {
#pragma omp barrier
#pragma omp master
    {
        VecSet(global_data.global_derivatives_at_node, 0.0);
        set_constant(global_data.derivatives_at_node, 0.0);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([&global_data](auto& elt) {
                auto& derivative = elt.data.derivative;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    subvector(global_data.derivatives_at_node,
                              derivative.local_nodeID[node] * GN::n_dddbath_terms,
                              GN::n_dddbath_terms) += derivative.dddbath_at_baryctr;
                }
            });
        }

        VecSetValuesBlocked(global_data.global_derivatives_at_node,
                            global_data.local_nodeIDs.size(),
                            &global_data.local_nodeIDs.front(),
                            global_data.derivatives_at_node.data(),
                            ADD_VALUES);

        VecAssemblyBegin(global_data.global_derivatives_at_node);
        VecAssemblyEnd(global_data.global_derivatives_at_node);

        VecScatterBegin(global_data.derivatives_scatter,
                        global_data.global_derivatives_at_node,
                        global_data.local_derivatives_at_node,
                        INSERT_VALUES,
                        SCATTER_FORWARD);
        VecScatterEnd(global_data.derivatives_scatter,
                      global_data.global_derivatives_at_node,
                      global_data.local_derivatives_at_node,
                      INSERT_VALUES,
                      SCATTER_FORWARD);

        double* d_ptr;
        VecGetArray(global_data.local_derivatives_at_node, &d_ptr);

        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            sim_units[su_id]->discretization.mesh.CallForEachElement([d_ptr](auto& elt) {
                auto& state      = elt.data.state[0];
                auto& derivative = elt.data.derivative;
                auto& internal   = elt.data.internal;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
                    for (uint dddbath = 0; dddbath < GN::n_dddbath_terms; ++dddbath) {
                        derivative.dddbath_lin(dddbath, node) =
                            d_ptr[GN::n_dddbath_terms * derivative.local_nodeID[node] + dddbath] /
                            derivative.node_mult[node];
                    }
                }
                state.dddbath          = elt.ProjectLinearToBasis(derivative.dddbath_lin);
                internal.dddbath_at_gp = elt.ComputeUgp(state.dddbath);
            });
        }

        VecRestoreArray(global_data.local_derivatives_at_node, &d_ptr);
    }
#pragma omp barrier
}
}
}

#endif
