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

    edge_state.q_hat = edge_bound.L2Projection(edge_internal.q_init_at_gp);

    /* Newton-Raphson iterator */

    uint iter = 0;
    while (true) {
        iter++;

        Problem::global_edge_boundary_iteration(stepper, edge_bound);

        if (iter == 100) {
            break;
        }

        double q_hat_norm     = norm(edge_state.q_hat);
        double delta_hat_norm = norm(edge_internal.rhs_global);

        if (delta_hat_norm / q_hat_norm < 1.0e-6) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
            edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
    }

    edge_bound.boundary.boundary_condition.ComputeNumericalFlux(edge_bound);
}

template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_iteration(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
            edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
    }

    /* Assemble global kernels */

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

    /* Increment q */

    for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat(SWE::Variables::ze, dof) += edge_internal.rhs_global[3 * dof];
        edge_state.q_hat(SWE::Variables::qx, dof) += edge_internal.rhs_global[3 * dof + 1];
        edge_state.q_hat(SWE::Variables::qy, dof) += edge_internal.rhs_global[3 * dof + 2];
    }
}
}
}

#endif