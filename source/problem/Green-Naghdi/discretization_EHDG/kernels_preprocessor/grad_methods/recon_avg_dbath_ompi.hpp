#ifndef RECON_AVG_DBATH_OMPI_HPP
#define RECON_AVG_DBATH_OMPI_HPP

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_dbath(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                       ProblemGlobalDataType& global_data,
                       const uint stage) {
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
                              GN::n_dimensions) += derivative.dbath_at_baryctr;
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
            sim_units[su_id]->discretization.mesh.CallForEachElement([stage, d_ptr](auto& elt) {
                auto& state      = elt.data.state[stage];
                auto& derivative = elt.data.derivative;
                auto& internal   = elt.data.internal;
                for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
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
    }
#pragma omp barrier
}

template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_ddbath(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                        ProblemGlobalDataType& global_data,
                        const uint stage) {
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
            sim_units[su_id]->discretization.mesh.CallForEachElement([stage, d_ptr](auto& elt) {
                auto& state      = elt.data.state[stage];
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
void reconstruct_dddbath(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                         ProblemGlobalDataType& global_data,
                         const uint stage) {
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
            sim_units[su_id]->discretization.mesh.CallForEachElement([stage, d_ptr](auto& elt) {
                auto& state      = elt.data.state[stage];
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
