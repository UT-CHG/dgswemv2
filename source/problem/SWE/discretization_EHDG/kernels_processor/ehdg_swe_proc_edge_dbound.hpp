#ifndef EHDG_SWE_PROC_EDGE_DBOUND_HPP
#define EHDG_SWE_PROC_EDGE_DBOUND_HPP

namespace SWE {
namespace EHDG {
template <typename EdgeDistributedType>
void Problem::global_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_dbound.boundary.boundary_condition.exchanger.GetEX(edge_dbound.boundary.boundary_condition.q_ex,
                                                            edge_dbound.boundary.boundary_condition.Fn_ex);

    edge_internal.q_init_at_gp = (boundary.q_at_gp + edge_dbound.boundary.boundary_condition.q_ex) / 2.0;

    edge_state.q_hat = edge_dbound.L2Projection(edge_internal.q_init_at_gp);

    /* Newton-Raphson iterator */

    uint iter = 0;
    while (true) {
        iter++;

        Problem::global_edge_distributed_iteration(stepper, edge_dbound);

        if (iter == 100) {
            break;
        }

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-8) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    edge_dbound.boundary.boundary_condition.ComputeNumericalFlux(edge_dbound);
}

template <typename EdgeDistributedType>
void Problem::global_edge_distributed_iteration(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    /* Assemble global kernels */

    edge_dbound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_dbound);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape_jacobian_vector<double, SWE::n_variables>(
                    edge_dbound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_dbound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

    /* Increment q */

    for (uint dof = 0; dof < edge_dbound.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat(SWE::Variables::ze, dof) += edge_internal.rhs_global[3 * dof];
        edge_state.q_hat(SWE::Variables::qx, dof) += edge_internal.rhs_global[3 * dof + 1];
        edge_state.q_hat(SWE::Variables::qy, dof) += edge_internal.rhs_global[3 * dof + 2];
    }
}
}
}

#endif