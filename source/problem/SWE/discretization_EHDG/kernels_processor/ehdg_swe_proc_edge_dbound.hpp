#ifndef EHDG_SWE_PROC_EDGE_DBOUND_HPP
#define EHDG_SWE_PROC_EDGE_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename EdgeDistributedType>
void Problem::global_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    /* Take average of in/ex state as initial trace state */
    Vector<double, SWE::n_variables> q_ex;
    Vector<double, SWE::n_variables> Fn_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(gp, q_ex, Fn_ex);

        edge_internal.q_avg_at_gp[gp] = (boundary.q_at_gp[gp] + q_ex) / 2.0;
    }

    edge_dbound.L2Projection(edge_internal.q_avg_at_gp, edge_state.q_hat);

    /* Newton-Raphson iterator */

    double q_norm      = 0;
    double q_norm_prev = 0;

    uint iter = 0;
    while (true) {
        Problem::global_edge_distributed_iteration(stepper, edge_dbound);

        iter++;

        if (iter == 100) {
            break;
        }

        q_norm      = 0;
        q_norm_prev = 0;

        for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
            q_norm += sqrNorm(edge_state.q_hat[dof]);
            q_norm_prev += sqrNorm(edge_state.q_hat_prev[dof]);
        }

        q_norm      = std::sqrt(q_norm);
        q_norm_prev = std::sqrt(q_norm_prev);

        if (std::abs(q_norm - q_norm_prev) / q_norm < 1.0e-6) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_dbound.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_internal.h_hat_at_gp[gp] = edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary.bath_at_gp[gp];
    }

    edge_dbound.boundary.boundary_condition.ComputeNumericalFlux(edge_dbound);
}

template <typename EdgeDistributedType>
void Problem::global_edge_distributed_iteration(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;
    auto& edge_global   = edge_dbound.edge_data.edge_global;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_dbound.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_internal.h_hat_at_gp[gp] = edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary.bath_at_gp[gp];
    }

    /* Assemble global kernels */

    edge_dbound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_dbound);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            Matrix<double, 3, 3> A =
                edge_dbound.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_hat_kernel_at_gp);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j)     = A(0, 0);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j) = A(1, 0);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j) = A(2, 0);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j + 1)     = A(0, 1);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j + 1) = A(1, 1);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j + 1) = A(2, 1);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j + 2)     = A(0, 2);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j + 2) = A(1, 2);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j + 2) = A(2, 2);
        }

        Vector<double, 3> b = edge_dbound.IntegrationLambda(dof_i, edge_global.rhs_kernel_at_gp);

        edge_global.global_rhs(3 * dof_i)     = b[0];
        edge_global.global_rhs(3 * dof_i + 1) = b[1];
        edge_global.global_rhs(3 * dof_i + 2) = b[2];
    }

    /* Solve global system for delta_q */

    edge_global.global_delta_q = edge_global.global_matrix.fullPivLu().solve(edge_global.global_rhs);

    /* Increment q */

    edge_state.q_hat_prev = edge_state.q_hat;

    for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat[dof][SWE::Variables::ze] += edge_global.global_delta_q(3 * dof);
        edge_state.q_hat[dof][SWE::Variables::qx] += edge_global.global_delta_q(3 * dof + 1);
        edge_state.q_hat[dof][SWE::Variables::qy] += edge_global.global_delta_q(3 * dof + 2);
    }
}
}
}

#endif