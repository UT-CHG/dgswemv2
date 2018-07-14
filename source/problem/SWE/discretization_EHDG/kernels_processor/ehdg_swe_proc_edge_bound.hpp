#ifndef EHDG_SWE_PROC_EDGE_BOUND_HPP
#define EHDG_SWE_PROC_EDGE_BOUND_HPP

namespace SWE {
namespace EHDG {
template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    /* Take in state as initial edge state */
    edge_internal.q_init_at_gp = boundary.q_at_gp;

    edge_bound.L2Projection(edge_internal.q_init_at_gp, edge_state.q_hat);

    /* Newton-Raphson iterator */

    double q_norm;
    double q_norm_prev;

    uint iter = 0;
    while (true) {
        Problem::global_edge_boundary_iteration(stepper, edge_bound);

        iter++;

        if (iter == 100) {
            break;
        }

        q_norm      = 0;
        q_norm_prev = 0;

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
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

    edge_bound.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    edge_bound.boundary.boundary_condition.ComputeNumericalFlux(edge_bound);
}

template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_iteration(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_bound.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    /* Assemble global kernels */

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
            blaze::submatrix(edge_internal.delta_hat_global,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) =
                edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp);
        }

        blaze::subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    int ipiv[edge_internal.rhs_global.size()];

    blaze::gesv(edge_internal.delta_hat_global, edge_internal.rhs_global, ipiv);

    /* Increment q */

    edge_state.q_hat_prev = edge_state.q_hat;

    for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat[dof][SWE::Variables::ze] += edge_internal.rhs_global[3 * dof];
        edge_state.q_hat[dof][SWE::Variables::qx] += edge_internal.rhs_global[3 * dof + 1];
        edge_state.q_hat[dof][SWE::Variables::qy] += edge_internal.rhs_global[3 * dof + 2];
    }
}
}
}

#endif