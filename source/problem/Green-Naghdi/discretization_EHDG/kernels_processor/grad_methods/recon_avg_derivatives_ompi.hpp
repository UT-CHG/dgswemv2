#ifndef RECON_AVG_DERIVATIVES_OMPI_HPP
#define RECON_AVG_DERIVATIVES_OMPI_HPP

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_dzedu(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, ProblemGlobalDataType& global_data) {
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
                              GN::n_dimensions) += derivative.dze_at_baryctr;
                    subvector(global_data.derivatives_at_node,
                              derivative.local_nodeID[node] * GN::n_dddbath_terms + GN::n_dimensions,
                              GN::n_du_terms) += derivative.du_at_baryctr;
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
                    for (uint dze = 0; dze < GN::n_dimensions; ++dze) {
                        derivative.dze_lin(dze, node) =
                            d_ptr[GN::n_dddbath_terms * derivative.local_nodeID[node] + dze] /
                            derivative.node_mult[node];
                    }
                    for (uint du = 0; du < GN::n_du_terms; ++du) {
                        derivative.du_lin(du, node) =
                            d_ptr[GN::n_dddbath_terms * derivative.local_nodeID[node] + GN::n_dimensions + du] /
                            derivative.node_mult[node];
                    }
                }
                state.dze         = elt.ProjectLinearToBasis(derivative.dze_lin);
                state.du          = elt.ProjectLinearToBasis(derivative.du_lin);
                internal.du_at_gp = elt.ComputeUgp(state.du);
            });
        }

        VecRestoreArray(global_data.local_derivatives_at_node, &d_ptr);
    }
#pragma omp barrier
}

template <typename OMPISimUnitType, typename ProblemGlobalDataType>
void reconstruct_ddu(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units, ProblemGlobalDataType& global_data) {
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
                              GN::n_ddu_terms) += derivative.ddu_at_baryctr;
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
                    for (uint ddu = 0; ddu < GN::n_ddu_terms; ++ddu) {
                        derivative.ddu_lin(ddu, node) =
                            d_ptr[GN::n_dddbath_terms * derivative.local_nodeID[node] + ddu] /
                            derivative.node_mult[node];
                    }
                }
                state.ddu          = elt.ProjectLinearToBasis(derivative.ddu_lin);
                internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
            });
        }

        VecRestoreArray(global_data.local_derivatives_at_node, &d_ptr);
    }
#pragma omp barrier
}
}
}

#endif
