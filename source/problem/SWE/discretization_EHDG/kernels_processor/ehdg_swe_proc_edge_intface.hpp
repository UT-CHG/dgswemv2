#ifndef EHDG_SWE_PROC_EDGE_INTFACE_HPP
#define EHDG_SWE_PROC_EDGE_INTFACE_HPP

namespace SWE {
namespace EHDG {
template <typename EdgeInterfaceType>
void Problem::global_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    /* Take average of in/ex state as initial trace state */
    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        edge_internal.q_avg_at_gp[gp] = (boundary_in.q_at_gp[gp] + boundary_ex.q_at_gp[gp_ex]) / 2.0;
    }

    edge_int.L2Projection(edge_internal.q_avg_at_gp, edge_state.q_hat);

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
            q_norm += sqrNorm(edge_state.q_hat[dof]);
            q_norm_prev += sqrNorm(edge_state.q_hat_prev[dof]);
        }

        q_norm      = std::sqrt(q_norm);
        q_norm_prev = std::sqrt(q_norm_prev);

        if ((std::abs(q_norm - q_norm_prev) / q_norm) < 1.0e-6) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_int.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_internal.h_hat_at_gp[gp] = edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary_in.bath_at_gp[gp];
    }

    edge_int.interface.specialization.ComputeNumericalFlux(edge_int);
}

template <typename EdgeInterfaceType>
void Problem::global_edge_interface_iteration(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;
    auto& edge_global   = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];

    edge_int.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_internal.h_hat_at_gp[gp] = edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary_in.bath_at_gp[gp];
    }

    /* Assemble global kernels */

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            Matrix<double, 3, 3> A = edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_hat_kernel_at_gp);

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

        Vector<double, 3> b = edge_int.IntegrationLambda(dof_i, edge_global.rhs_kernel_at_gp);

        edge_global.global_rhs(3 * dof_i)     = b[0];
        edge_global.global_rhs(3 * dof_i + 1) = b[1];
        edge_global.global_rhs(3 * dof_i + 2) = b[2];
    }

    /* Solve global system for delta_q */

    edge_global.global_delta_q = edge_global.global_matrix.fullPivLu().solve(edge_global.global_rhs);

    /* Increment q */

    edge_state.q_hat_prev = edge_state.q_hat;

    for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat[dof][SWE::Variables::ze] += edge_global.global_delta_q(3 * dof);
        edge_state.q_hat[dof][SWE::Variables::qx] += edge_global.global_delta_q(3 * dof + 1);
        edge_state.q_hat[dof][SWE::Variables::qy] += edge_global.global_delta_q(3 * dof + 2);
    }
}
}
}

#endif