#ifndef IHDG_SWE_PROC_EDGE_INTFACE_HPP
#define IHDG_SWE_PROC_EDGE_INTFACE_HPP

#include <eigen3/Eigen/Dense>

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
void Problem::local_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& internal_in = edge_int.interface.data_in.internal;
    auto& internal_ex = edge_int.interface.data_ex.internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);

    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_intface_LF(edge_int);

    // Add (dtau/dq_hat * del_q - tau * del_q) term to dF_hat_dq_hat
    // and tau term to dF_hat_dq
    add_dF_hat_tau_terms_intface_LF(edge_int);

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); dof_j++) {
            submatrix(internal_in.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.interface.IntegrationPhiPhiIN(dof_i, dof_j, boundary_in.dF_hat_dq_at_gp));
        }

        subvector(internal_in.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -edge_int.interface.IntegrationPhiIN(dof_i, boundary_in.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); dof_j++) {
            submatrix(internal_ex.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.interface.IntegrationPhiPhiEX(dof_i, dof_j, boundary_ex.dF_hat_dq_at_gp));
        }

        subvector(internal_ex.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -edge_int.interface.IntegrationPhiEX(dof_i, boundary_ex.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); dof_j++) {
            submatrix(boundary_in.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.IntegrationPhiLambdaIN(dof_i, dof_j, boundary_in.dF_hat_dq_hat_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); dof_j++) {
            submatrix(boundary_ex.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.IntegrationPhiLambdaEX(dof_i, dof_j, boundary_ex.dF_hat_dq_hat_at_gp));
        }
    }
}

template <typename EdgeInterfaceType>
void Problem::global_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); dof_j++) {
            submatrix(boundary_in.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.IntegrationPhiLambdaIN(dof_j, dof_i, boundary_in.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); dof_j++) {
            submatrix(boundary_ex.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.IntegrationPhiLambdaEX(dof_j, dof_i, boundary_ex.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); dof_j++) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_int.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }
}
}
}

#endif