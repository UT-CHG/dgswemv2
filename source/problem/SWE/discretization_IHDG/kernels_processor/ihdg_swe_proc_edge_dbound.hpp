#ifndef IHDG_SWE_PROC_EDGE_DBOUND_HPP
#define IHDG_SWE_PROC_EDGE_DBOUND_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename EdgeDistributedType>
void Problem::init_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound) {
    if (stepper.GetOrder() == 2) {
        auto& edge_state    = edge_dbound.edge_data.edge_state;
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& internal = edge_dbound.boundary.data.internal;
        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        // Add tau * del_q terms to F_hat
        add_F_hat_tau_terms_bound_LF(edge_dbound);

        double theta = stepper.GetTheta();

        for (uint dof_i = 0; dof_i < edge_dbound.boundary.data.get_ndof(); ++dof_i) {
            subvector(internal.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) +=
                -theta * edge_dbound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp);
        }
    }
}

template <typename StepperType, typename EdgeDistributedType>
void Problem::local_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& internal = edge_dbound.boundary.data.internal;
    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_bound_LF(edge_dbound);

    // Add (dtau/dq_hat * del_q - tau) term to dF_hat_dq_hat
    // and tau term to dF_hat_dq
    add_dF_hat_tau_terms_bound_LF(edge_dbound);

    double theta = stepper.GetTheta();

    for (uint dof_i = 0; dof_i < edge_dbound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_dbound.boundary.IntegrationPhiPhi(dof_j, dof_i, boundary.dF_hat_dq_at_gp));
        }

        subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -(1.0 - theta) * edge_dbound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_dbound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            submatrix(boundary.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_dbound.IntegrationPhiLambda(dof_i, dof_j, boundary.dF_hat_dq_hat_at_gp));
        }
    }
}

template <typename StepperType, typename EdgeDistributedType>
void Problem::global_edge_distributed_kernel(const StepperType& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_dbound.boundary.boundary_condition.ComputeGlobalKernels(edge_dbound);

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(boundary.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_dbound.IntegrationPhiLambda(dof_j, dof_i, boundary.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_dbound.IntegrationLambdaLambda(dof_j, dof_i, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_dbound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }
}
}
}

#endif