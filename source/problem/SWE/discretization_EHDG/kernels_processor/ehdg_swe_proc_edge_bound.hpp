#ifndef EHDG_SWE_PROC_EDGE_BOUND_HPP
#define EHDG_SWE_PROC_EDGE_BOUND_HPP

namespace SWE {
namespace EHDG {
template <typename StepperType, typename EdgeBoundaryType>
void Problem::global_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    /* Newton-Raphson iterator */

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        Problem::global_edge_boundary_iteration(stepper, edge_bound);

        double delta_hat_norm = norm(edge_bound.edge_data.edge_internal.rhs_global);

        if (delta_hat_norm < 1.0e-12) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_bound.boundary.boundary_condition.ComputeNumericalFlux(edge_bound);
}

template <typename StepperType, typename EdgeBoundaryType>
void Problem::global_edge_boundary_iteration(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

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
                reshape<double, SWE::n_variables>(
                    edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

    /* Increment q */

    edge_state.q_hat +=
        reshape<double, SWE::n_variables, SO::ColumnMajor>(edge_internal.rhs_global, edge_bound.edge_data.get_ndof());
}
}
}

#endif