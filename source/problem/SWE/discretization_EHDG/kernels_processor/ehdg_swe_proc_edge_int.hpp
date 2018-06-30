#ifndef EHDG_SWE_PROC_EDGE_INT_HPP
#define EHDG_SWE_PROC_EDGE_INT_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_EHDG/stabilization_parameters/ehdg_swe_stabilization_params.hpp"

namespace SWE {
namespace EHDG {
template <typename EdgeInterfaceType>
void Problem::global_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    /* Take average of in/ex state as initial trace state */
    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        edge_state.ze_avg_at_gp[gp] = (boundary_in.ze_at_gp[gp] + boundary_ex.ze_at_gp[gp_ex]) / 2.0;
        edge_state.qx_avg_at_gp[gp] = (boundary_in.qx_at_gp[gp] + boundary_ex.qx_at_gp[gp_ex]) / 2.0;
        edge_state.qy_avg_at_gp[gp] = (boundary_in.qy_at_gp[gp] + boundary_ex.qy_at_gp[gp_ex]) / 2.0;
    }

    edge_int.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_int.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_int.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);

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
            q_norm += std::pow(edge_state.ze_hat[dof], 2) + std::pow(edge_state.qx_hat[dof], 2) +
                      std::pow(edge_state.qy_hat[dof], 2);

            q_norm_prev += std::pow(edge_state.ze_hat_prev[dof], 2) + std::pow(edge_state.qx_hat_prev[dof], 2) +
                           std::pow(edge_state.qy_hat_prev[dof], 2);
        }

        q_norm      = std::sqrt(q_norm);
        q_norm_prev = std::sqrt(q_norm_prev);

        if ((std::abs(q_norm - q_norm_prev) / q_norm) < 1.0e-6) {
            break;
        }
    }

    Problem::edge_interface_compute_flux(stepper, edge_int);
}

template <typename EdgeInterfaceType>
void Problem::global_edge_interface_iteration(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary_in.bath_at_gp[gp];
    }

    /* Assemble kernels for global matrix */

    // Add tau terms
    add_delta_kernel_tau_terms_LF(edge_int);

    /* Assemble kernels for global rhs vector */
    double nx_in, ny_in;
    double nx_ex, ny_ex;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        // Add F_in * n_in terms
        nx_in = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny_in = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

        edge_global.ze_rhs_kernel_at_gp[gp] = -(boundary_in.ze_flux_at_gp[GlobalCoord::x][gp] * nx_in +
                                                boundary_in.ze_flux_at_gp[GlobalCoord::y][gp] * ny_in);
        edge_global.qx_rhs_kernel_at_gp[gp] = -(boundary_in.qx_flux_at_gp[GlobalCoord::x][gp] * nx_in +
                                                boundary_in.qx_flux_at_gp[GlobalCoord::y][gp] * ny_in);
        edge_global.qy_rhs_kernel_at_gp[gp] = -(boundary_in.qy_flux_at_gp[GlobalCoord::x][gp] * nx_in +
                                                boundary_in.qy_flux_at_gp[GlobalCoord::y][gp] * ny_in);

        // Add F_ex * n_ex terms
        nx_ex = edge_int.interface.surface_normal_ex[gp_ex][GlobalCoord::x];
        ny_ex = edge_int.interface.surface_normal_ex[gp_ex][GlobalCoord::y];

        edge_global.ze_rhs_kernel_at_gp[gp] += -(boundary_ex.ze_flux_at_gp[GlobalCoord::x][gp_ex] * nx_ex +
                                                 boundary_ex.ze_flux_at_gp[GlobalCoord::y][gp_ex] * ny_ex);
        edge_global.qx_rhs_kernel_at_gp[gp] += -(boundary_ex.qx_flux_at_gp[GlobalCoord::x][gp_ex] * nx_ex +
                                                 boundary_ex.qx_flux_at_gp[GlobalCoord::y][gp_ex] * ny_ex);
        edge_global.qy_rhs_kernel_at_gp[gp] += -(boundary_ex.qy_flux_at_gp[GlobalCoord::x][gp_ex] * nx_ex +
                                                 boundary_ex.qy_flux_at_gp[GlobalCoord::y][gp_ex] * ny_ex);
    }

    // Add tau terms
    add_rhs_kernel_tau_terms_LF(edge_int);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            edge_global.global_matrix(3 * dof_i, 3 * dof_j) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_ze_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_ze_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_ze_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j + 1) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qx_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j + 1) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qx_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j + 1) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qx_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_i, 3 * dof_j + 2) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qy_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_i + 1, 3 * dof_j + 2) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qy_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_i + 2, 3 * dof_j + 2) =
                edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_global.delta_qy_hat_kernel_at_gp[Variables::qy]);
        }

        edge_global.global_rhs(3 * dof_i)     = edge_int.IntegrationLambda(dof_i, edge_global.ze_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_i + 1) = edge_int.IntegrationLambda(dof_i, edge_global.qx_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_i + 2) = edge_int.IntegrationLambda(dof_i, edge_global.qy_rhs_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    edge_global.global_delta_q = edge_global.global_matrix.fullPivLu().solve(edge_global.global_rhs);

    /* Increment q */

    for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
        edge_state.ze_hat_prev[dof] = edge_state.ze_hat[dof];
        edge_state.qx_hat_prev[dof] = edge_state.qx_hat[dof];
        edge_state.qy_hat_prev[dof] = edge_state.qy_hat[dof];

        edge_state.ze_hat[dof] += edge_global.global_delta_q(3 * dof);
        edge_state.qx_hat[dof] += edge_global.global_delta_q(3 * dof + 1);
        edge_state.qy_hat[dof] += edge_global.global_delta_q(3 * dof + 2);
    }
}

template <typename EdgeInterfaceType>
void Problem::edge_interface_compute_flux(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_int.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary_in.bath_at_gp[gp];
    }

    double nx_in, ny_in;
    double nx_ex, ny_ex;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        nx_in = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny_in = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

        boundary_in.ze_numerical_flux_at_gp[gp] = boundary_in.ze_flux_at_gp[GlobalCoord::x][gp] * nx_in +
                                                  boundary_in.ze_flux_at_gp[GlobalCoord::y][gp] * ny_in;
        boundary_in.qx_numerical_flux_at_gp[gp] = boundary_in.qx_flux_at_gp[GlobalCoord::x][gp] * nx_in +
                                                  boundary_in.qx_flux_at_gp[GlobalCoord::y][gp] * ny_in;
        boundary_in.qy_numerical_flux_at_gp[gp] = boundary_in.qy_flux_at_gp[GlobalCoord::x][gp] * nx_in +
                                                  boundary_in.qy_flux_at_gp[GlobalCoord::y][gp] * ny_in;

        nx_ex = edge_int.interface.surface_normal_ex[gp_ex][GlobalCoord::x];
        ny_ex = edge_int.interface.surface_normal_ex[gp_ex][GlobalCoord::y];

        boundary_ex.ze_numerical_flux_at_gp[gp_ex] = boundary_ex.ze_flux_at_gp[GlobalCoord::x][gp_ex] * nx_ex +
                                                     boundary_ex.ze_flux_at_gp[GlobalCoord::y][gp_ex] * ny_ex;
        boundary_ex.qx_numerical_flux_at_gp[gp_ex] = boundary_ex.qx_flux_at_gp[GlobalCoord::x][gp_ex] * nx_ex +
                                                     boundary_ex.qx_flux_at_gp[GlobalCoord::y][gp_ex] * ny_ex;
        boundary_ex.qy_numerical_flux_at_gp[gp_ex] = boundary_ex.qy_flux_at_gp[GlobalCoord::x][gp_ex] * nx_ex +
                                                     boundary_ex.qy_flux_at_gp[GlobalCoord::y][gp_ex] * ny_ex;
    }

    // Add tau terms
    add_flux_tau_terms_intface_LF(edge_int);
}
}
}

#endif