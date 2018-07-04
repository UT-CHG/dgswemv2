#ifndef IHDG_SWE_PROC_EDGE_DBOUND_HPP
#define IHDG_SWE_PROC_EDGE_DBOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeDistributedType>
void Problem::global_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state = edge_dbound.edge_data.edge_state;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    /* Take average of in/ex state as initial trace state */
    double ze_ex, qx_ex, qy_ex;
    double ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(
            gp, ze_ex, qx_ex, qy_ex, ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex);

        edge_state.ze_avg_at_gp[gp] = (boundary.ze_at_gp[gp] + ze_ex) / 2.0;
        edge_state.qx_avg_at_gp[gp] = (boundary.qx_at_gp[gp] + qx_ex) / 2.0;
        edge_state.qy_avg_at_gp[gp] = (boundary.qy_at_gp[gp] + qy_ex) / 2.0;
    }

    edge_dbound.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_dbound.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_dbound.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);

    /* Newton-Raphson iterator */

    double q_norm      = 0;
    double q_norm_prev = 0;

    uint iter = 0;
    while (true) {
        Problem::global_edge_boundary_iteration(stepper, edge_dbound);

        iter++;

        if (iter == 100) {
            break;
        }

        q_norm      = 0;
        q_norm_prev = 0;

        for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
            q_norm += std::pow(edge_state.ze_hat[dof], 2) + std::pow(edge_state.qx_hat[dof], 2) +
                      std::pow(edge_state.qy_hat[dof], 2);

            q_norm_prev += std::pow(edge_state.ze_hat_prev[dof], 2) + std::pow(edge_state.qx_hat_prev[dof], 2) +
                           std::pow(edge_state.qy_hat_prev[dof], 2);
        }

        q_norm      = std::sqrt(q_norm);
        q_norm_prev = std::sqrt(q_norm_prev);

        if (std::abs(q_norm - q_norm_prev) / q_norm < 1.0e-6) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_dbound.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_dbound.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_dbound.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary.bath_at_gp[gp];
    }

    edge_dbound.boundary.boundary_condition.ComputeNumericalFlux(edge_dbound);
}

template <typename EdgeDistributedType>
void Problem::global_edge_distributed_iteration(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state  = edge_dbound.edge_data.edge_state;
    auto& edge_global = edge_dbound.edge_data.edge_global;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_dbound.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_dbound.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_dbound.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary.bath_at_gp[gp];
    }

    /* Assemble global kernels */

    edge_dbound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_dbound);

    /* Assemble global system */

    for (uint dof_row = 0; dof_row < edge_dbound.edge_data.get_ndof(); ++dof_row) {
        for (uint dof_col = 0; dof_col < edge_dbound.edge_data.get_ndof(); ++dof_col) {
            edge_global.global_matrix(3 * dof_row, 3 * dof_col) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_ze_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_row + 1, 3 * dof_col) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_ze_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_row + 2, 3 * dof_col) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_ze_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_row, 3 * dof_col + 1) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qx_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_row + 1, 3 * dof_col + 1) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qx_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_row + 2, 3 * dof_col + 1) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qx_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_row, 3 * dof_col + 2) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qy_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_row + 1, 3 * dof_col + 2) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qy_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_row + 2, 3 * dof_col + 2) = edge_dbound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qy_hat_kernel_at_gp[Variables::qy]);
        }

        edge_global.global_rhs(3 * dof_row) = edge_dbound.IntegrationLambda(dof_row, edge_global.ze_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_row + 1) =
            edge_dbound.IntegrationLambda(dof_row, edge_global.qx_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_row + 2) =
            edge_dbound.IntegrationLambda(dof_row, edge_global.qy_rhs_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    edge_global.global_delta_q = edge_global.global_matrix.fullPivLu().solve(edge_global.global_rhs);

    /* Increment q */

    for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
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