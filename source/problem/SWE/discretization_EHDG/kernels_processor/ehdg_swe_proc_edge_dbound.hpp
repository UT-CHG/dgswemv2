#ifndef EHDG_SWE_PROC_EDGE_DBOUND_HPP
#define EHDG_SWE_PROC_EDGE_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename StepperType, typename EdgeDistributedType>
void Problem::global_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    /* Newton-Raphson iterator */

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        Problem::global_edge_distributed_iteration(stepper, edge_dbound);

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-12) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    edge_dbound.boundary.boundary_condition.ComputeNumericalFlux(edge_dbound);
}

template <typename StepperType, typename EdgeDistributedType>
void Problem::global_edge_distributed_iteration(const StepperType& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    /* Assemble global kernels */

    edge_dbound.boundary.boundary_condition.ComputeGlobalKernels(edge_dbound);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_dbound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_dbound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

    /* Increment q */

    edge_state.q_hat +=
        reshape<double, SWE::n_variables, SO::ColumnMajor>(edge_internal.rhs_global, edge_dbound.edge_data.get_ndof());
}
}
}

#endif