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
    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        edge_internal.q_init_at_gp[gp] = (boundary_in.q_at_gp[gp] + boundary_ex.q_at_gp[gp_ex]) / 2.0;
    }

    edge_int.L2Projection(edge_internal.q_init_at_gp, edge_state.q_hat);

    /* Newton-Raphson iterator */

    uint iter = 0;
    while (true) {
        iter++;

        Problem::global_edge_interface_iteration(stepper, edge_int);

        if (iter == 100) {
            break;
        }

        double q_hat_norm     = 0;
        double delta_hat_norm = blaze::norm(edge_internal.rhs_global);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            q_hat_norm += sqrNorm(edge_state.q_hat[dof]);
        }

        q_hat_norm = std::sqrt(q_hat_norm);

        if (delta_hat_norm / q_hat_norm < 1.0e-6) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_int.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary_in.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    edge_int.interface.specialization.ComputeNumericalFlux(edge_int);
}

template <typename EdgeInterfaceType>
void Problem::global_edge_interface_iteration(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];

    edge_int.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary_in.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    /* Assemble global kernels */

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            blaze::submatrix(edge_internal.delta_hat_global,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp);
        }

        blaze::subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_int.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    int ipiv[edge_internal.rhs_global.size()];

    blaze::gesv(edge_internal.delta_hat_global, edge_internal.rhs_global, ipiv);

    /* Increment q */

    edge_state.q_hat_prev = edge_state.q_hat;

    for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat[dof][SWE::Variables::ze] += edge_internal.rhs_global[3 * dof];
        edge_state.q_hat[dof][SWE::Variables::qx] += edge_internal.rhs_global[3 * dof + 1];
        edge_state.q_hat[dof][SWE::Variables::qy] += edge_internal.rhs_global[3 * dof + 2];
    }
}
}
}

#endif