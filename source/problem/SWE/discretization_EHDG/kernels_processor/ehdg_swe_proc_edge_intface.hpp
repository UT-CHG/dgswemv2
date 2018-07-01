#ifndef EHDG_SWE_PROC_EDGE_INTFACE_HPP
#define EHDG_SWE_PROC_EDGE_INTFACE_HPP

#include <eigen3/Eigen/Dense>

namespace SWE {
namespace EHDG {
template <typename EdgeInterfaceType>
void Problem::global_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    /* Take average of in/ex state as initial trace state */
    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        edge_state.ze_avg_at_gp[gp] = (boundary_in.ze_at_gp[gp] + boundary_ex.ze_at_gp[gp_ex]) / 2.0;
        edge_state.qx_avg_at_gp[gp] = (boundary_in.qx_at_gp[gp] + boundary_ex.qx_at_gp[gp_ex]) / 2.0;
        edge_state.qy_avg_at_gp[gp] = (boundary_in.qy_at_gp[gp] + boundary_ex.qy_at_gp[gp_ex]) / 2.0;
    }

    edge_int.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_int.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_int.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);

    /* Newton-Raphson iterator */

    double q_norm;
    double q_norm_prev;

    uint iter = 0;
    while (true) {
        Problem::global_edge_interface_iteration(stepper, edge_int);

        iter++;

        if (iter == 100) {
            break;
        }

        q_norm      = 0;
        q_norm_prev = 0;

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            q_norm += std::pow(edge_state.ze_hat[dof], 2) + std::pow(edge_state.qx_hat[dof], 2) +
                      std::pow(edge_state.qy_hat[dof], 2);

            q_norm_prev += std::pow(edge_state.ze_hat_prev[dof], 2) + std::pow(edge_state.qx_hat_prev[dof], 2) +
                           std::pow(edge_state.qy_hat_prev[dof], 2);
        }

        q_norm      = std::sqrt(q_norm);
        q_norm_prev = std::sqrt(q_norm_prev);

        if ((std::abs(q_norm - q_norm_prev) / q_norm) < 1.0e-6) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_int.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary_in.bath_at_gp[gp];
    }

    edge_int.interface.specialization.ComputeNumericalFlux(edge_int);
}

template <typename EdgeInterfaceType>
void Problem::global_edge_interface_iteration(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary_in.bath_at_gp[gp];
    }

    /* Assemble global kernels */

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            edge_global.global_matrix(3 * dof_i, 3 * dof_j) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_ze_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_ze_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_ze_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j + 1) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qx_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j + 1) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qx_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j + 1) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qx_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j + 2) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qy_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j + 2) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qy_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j + 2) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qy_hat_kernel_at_gp[Variables::qy]);
        }

        edge_global.global_rhs(3 * dof_i)     = edge_int.IntegrationLambda(dof_i, edge_global.ze_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_i + 1) = edge_int.IntegrationLambda(dof_i, edge_global.qx_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_i + 2) = edge_int.IntegrationLambda(dof_i, edge_global.qy_rhs_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    edge_global.global_delta_q = edge_global.global_matrix.fullPivLu().solve(edge_global.global_rhs);

    /* Increment q */

    for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
        edge_state.ze_hat_prev[dof] = edge_state.ze_hat[dof];
        edge_state.qx_hat_prev[dof] = edge_state.qx_hat[dof];
        edge_state.qy_hat_prev[dof] = edge_state.qy_hat[dof];

        edge_state.ze_hat[dof] += edge_global.global_delta_q(3 * dof);
        edge_state.qx_hat[dof] += edge_global.global_delta_q(3 * dof + 1);
        edge_state.qy_hat[dof] += edge_global.global_delta_q(3 * dof + 2);
    }
}
}
}

#endif