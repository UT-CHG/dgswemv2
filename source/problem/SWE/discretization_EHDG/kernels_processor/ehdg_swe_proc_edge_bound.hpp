#ifndef EHDG_SWE_PROC_EDGE_BOUND_HPP
#define EHDG_SWE_PROC_EDGE_BOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_EHDG/stabilization_parameters/ehdg_swe_stabilization_params.hpp"

namespace SWE {
namespace EHDG {
template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    const uint stage = stepper.GetStage();

    auto& edge_state  = edge_bound.edge_data.edge_state;
    auto& edge_global = edge_bound.edge_data.edge_global;

    auto& state    = edge_bound.data.state[stage];
    auto& boundary = edge_bound.data.boundary[edge_bound.bound_id];

    /* Take linear combination of in/ex state as initial edge state */
    /* THIS IS PERIODIC BC FOR MANU SOLUTION */
    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_state.ze_avg_at_gp[gp] = boundary.ze_at_gp[gp];
        edge_state.qx_avg_at_gp[gp] = boundary.qx_at_gp[gp];
        edge_state.qy_avg_at_gp[gp] = boundary.qy_at_gp[gp];
    }
    /* THIS IS PERIODIC BC FOR MANU SOLUTION */

    edge_bound.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_bound.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_bound.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);

    /*double q_norm      = 0;
    double q_norm_prev = 0;

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
            q_norm += std::pow(edge_state.ze_hat[dof], 2) + std::pow(edge_state.qx_hat[dof], 2) +
                      std::pow(edge_state.qy_hat[dof], 2);

            q_norm_prev += std::pow(edge_state.ze_hat_prev[dof], 2) + std::pow(edge_state.qx_hat_prev[dof], 2) +
                           std::pow(edge_state.qy_hat_prev[dof], 2);
        }

        q_norm      = std::sqrt(q_norm);
        q_norm_prev = std::sqrt(q_norm_prev);

        if (std::abs(q_norm - q_norm_prev) / q_norm < 1.0e-3) {
            break;
        }
    }*/

    Problem::edge_boundary_compute_flux(stepper, edge_bound);
}

template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_iteration(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    const uint stage = stepper.GetStage();

    auto& edge_state  = edge_bound.edge_data.edge_state;
    auto& edge_global = edge_bound.edge_data.edge_global;

    auto& state    = edge_bound.data.state[stage];
    auto& boundary = edge_bound.data.boundary[edge_bound.bound_id];

    edge_bound.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_bound.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_bound.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    /* THIS IS PERIODIC BC FOR MANU SOLUTION */
    std::vector<double> ze_at_gp_ex(edge_bound.edge_data.get_ngp());
    std::vector<double> qx_at_gp_ex(edge_bound.edge_data.get_ngp());
    std::vector<double> qy_at_gp_ex(edge_bound.edge_data.get_ngp());

    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_bound.edge_data.get_ngp() - gp - 1;

        ze_at_gp_ex[gp_ex] = boundary.ze_at_gp[gp];
        qx_at_gp_ex[gp_ex] = boundary.qx_at_gp[gp];
        qy_at_gp_ex[gp_ex] = boundary.qy_at_gp[gp];
    }
    /* THIS IS PERIODIC BC FOR MANU SOLUTION */

    /* Assemble kernels for global matrix */

    double nx = 0.0;
    double ny = 0.0;

    double u_hat = 0.0;
    double v_hat = 0.0;

    double uuh_hat = 0.0;
    double vvh_hat = 0.0;
    double uvh_hat = 0.0;
    double pe_hat  = 0.0;

    gp_ex = 0;
    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_bound.edge_data.get_ngp() - gp - 1;

        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary.bath_at_gp[gp];

        nx = edge_bound.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.surface_normal[gp][GlobalCoord::y];

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        // Start with dF_hat/dq_hat * n_in terms
        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] = 0.0;
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qx][gp] =
            (-u_hat * u_hat + Global::g * edge_state.h_hat_at_gp[gp]) * nx - u_hat * v_hat * ny;
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qy][gp] =
            -u_hat * v_hat * nx + (-v_hat * v_hat + Global::g * edge_state.h_hat_at_gp[gp]) * ny;

        edge_global.delta_qx_hat_kernel_at_gp[Variables::ze][gp] = nx;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] = 2 * u_hat * nx + v_hat * ny;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qy][gp] = v_hat * nx;

        edge_global.delta_qy_hat_kernel_at_gp[Variables::ze][gp] = ny;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qx][gp] = u_hat * ny;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] = u_hat * nx + 2 * v_hat * ny;

        // Add dF_hat/dq_hat * n_ex terms
        nx = -edge_bound.surface_normal[gp][GlobalCoord::x];
        ny = -edge_bound.surface_normal[gp][GlobalCoord::y];

        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] += 0.0;
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qx][gp] +=
            (-u_hat * u_hat + Global::g * edge_state.h_hat_at_gp[gp]) * nx - u_hat * v_hat * ny;
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qy][gp] +=
            -u_hat * v_hat * nx + (-v_hat * v_hat + Global::g * edge_state.h_hat_at_gp[gp]) * ny;

        edge_global.delta_qx_hat_kernel_at_gp[Variables::ze][gp] += nx;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] += 2 * u_hat * nx + v_hat * ny;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qy][gp] += v_hat * nx;

        edge_global.delta_qy_hat_kernel_at_gp[Variables::ze][gp] += ny;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qx][gp] += u_hat * ny;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] += u_hat * nx + 2 * v_hat * ny;
    }

    // Add tau terms
    add_delta_kernel_tau_LF(edge_bound.surface_normal,
                            edge_state.ze_hat_at_gp,
                            edge_state.qx_hat_at_gp,
                            edge_state.qy_hat_at_gp,
                            edge_state.h_hat_at_gp,
                            edge_global.delta_ze_hat_kernel_at_gp,
                            edge_global.delta_qx_hat_kernel_at_gp,
                            edge_global.delta_qy_hat_kernel_at_gp);

    // Add dtau/dq_hat * diff_q terms
    add_delta_kernel_partial_tau_diff_q_LF(edge_bound.surface_normal,
                                           boundary.ze_at_gp,
                                           boundary.qx_at_gp,
                                           boundary.qy_at_gp,
                                           ze_at_gp_ex,
                                           qx_at_gp_ex,
                                           qy_at_gp_ex,
                                           edge_state.ze_hat_at_gp,
                                           edge_state.qx_hat_at_gp,
                                           edge_state.qy_hat_at_gp,
                                           edge_state.h_hat_at_gp,
                                           edge_global.delta_ze_hat_kernel_at_gp,
                                           edge_global.delta_qx_hat_kernel_at_gp,
                                           edge_global.delta_qy_hat_kernel_at_gp);

    /* Assemble kernels for global rhs vector */

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_bound.edge_data.get_ngp() - gp - 1;

        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary.bath_at_gp[gp];

        nx = edge_bound.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.surface_normal[gp][GlobalCoord::y];

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        uuh_hat = u_hat * edge_state.qx_hat_at_gp[gp];
        vvh_hat = v_hat * edge_state.qy_hat_at_gp[gp];
        uvh_hat = u_hat * edge_state.qy_hat_at_gp[gp];
        pe_hat  = Global::g * (0.5 * std::pow(edge_state.ze_hat_at_gp[gp], 2) +
                              edge_state.ze_hat_at_gp[gp] * boundary.bath_at_gp[gp]);

        // Add F_hat * n_in terms
        edge_global.ze_rhs_kernel_at_gp[gp] = -(edge_state.qx_hat_at_gp[gp] * nx + edge_state.qy_hat_at_gp[gp] * ny);
        edge_global.qx_rhs_kernel_at_gp[gp] = -((uuh_hat + pe_hat) * nx + uvh_hat * ny);
        edge_global.qy_rhs_kernel_at_gp[gp] = -(uvh_hat * nx + (vvh_hat + pe_hat) * ny);

        // Add F_hat * n_ex terms
        nx = -edge_bound.surface_normal[gp][GlobalCoord::x];
        ny = -edge_bound.surface_normal[gp][GlobalCoord::y];

        edge_global.ze_rhs_kernel_at_gp[gp] -= edge_state.qx_hat_at_gp[gp] * nx + edge_state.qy_hat_at_gp[gp] * ny;
        edge_global.qx_rhs_kernel_at_gp[gp] -= (uuh_hat + pe_hat) * nx + uvh_hat * ny;
        edge_global.qy_rhs_kernel_at_gp[gp] -= uvh_hat * nx + (vvh_hat + pe_hat) * ny;
    }

    // Add tau * diff_q terms
    add_rhs_kernel_tau_diff_q_LF(edge_bound.surface_normal,
                                 boundary.ze_at_gp,
                                 boundary.qx_at_gp,
                                 boundary.qy_at_gp,
                                 ze_at_gp_ex,
                                 qx_at_gp_ex,
                                 qy_at_gp_ex,
                                 edge_state.ze_hat_at_gp,
                                 edge_state.qx_hat_at_gp,
                                 edge_state.qy_hat_at_gp,
                                 edge_state.h_hat_at_gp,
                                 edge_global.ze_rhs_kernel_at_gp,
                                 edge_global.qx_rhs_kernel_at_gp,
                                 edge_global.qy_rhs_kernel_at_gp);

    /* Assemble global system */

    for (uint dof_row = 0; dof_row < edge_bound.edge_data.get_ndof(); ++dof_row) {
        for (uint dof_col = 0; dof_col < edge_bound.edge_data.get_ndof(); ++dof_col) {
            edge_global.global_matrix(3 * dof_row, 3 * dof_col) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_ze_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_row + 1, 3 * dof_col) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_ze_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_row + 2, 3 * dof_col) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_ze_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_row, 3 * dof_col + 1) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qx_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_row + 1, 3 * dof_col + 1) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qx_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_row + 2, 3 * dof_col + 1) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qx_hat_kernel_at_gp[Variables::qy]);

            edge_global.global_matrix(3 * dof_row, 3 * dof_col + 2) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qy_hat_kernel_at_gp[Variables::ze]);
            edge_global.global_matrix(3 * dof_row + 1, 3 * dof_col + 2) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qy_hat_kernel_at_gp[Variables::qx]);
            edge_global.global_matrix(3 * dof_row + 2, 3 * dof_col + 2) = edge_bound.IntegrationLambdaLambda(
                dof_row, dof_col, edge_global.delta_qy_hat_kernel_at_gp[Variables::qy]);
        }

        edge_global.global_rhs(3 * dof_row) = edge_bound.IntegrationLambda(dof_row, edge_global.ze_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_row + 1) =
            edge_bound.IntegrationLambda(dof_row, edge_global.qx_rhs_kernel_at_gp);
        edge_global.global_rhs(3 * dof_row + 2) =
            edge_bound.IntegrationLambda(dof_row, edge_global.qy_rhs_kernel_at_gp);
    }

    /* Solve global systel for delta_q */

    edge_global.global_delta_q = edge_global.global_matrix.fullPivLu().solve(edge_global.global_rhs);

    /* Increment q */

    for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
        edge_state.ze_hat_prev[dof] = edge_state.ze_hat[dof];
        edge_state.qx_hat_prev[dof] = edge_state.qx_hat[dof];
        edge_state.qy_hat_prev[dof] = edge_state.qy_hat[dof];

        edge_state.ze_hat[dof] += edge_global.global_delta_q(3 * dof);
        edge_state.qx_hat[dof] += edge_global.global_delta_q(3 * dof + 1);
        edge_state.qy_hat[dof] += edge_global.global_delta_q(3 * dof + 2);
    }
}

template <typename EdgeBoundaryType>
void Problem::edge_boundary_compute_flux(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    const uint stage = stepper.GetStage();

    auto& edge_state  = edge_bound.edge_data.edge_state;
    auto& edge_global = edge_bound.edge_data.edge_global;

    auto& state    = edge_bound.data.state[stage];
    auto& boundary = edge_bound.data.boundary[edge_bound.bound_id];

    edge_bound.ComputeUgp(edge_state.ze_hat, edge_state.ze_hat_at_gp);
    edge_bound.ComputeUgp(edge_state.qx_hat, edge_state.qx_hat_at_gp);
    edge_bound.ComputeUgp(edge_state.qy_hat, edge_state.qy_hat_at_gp);

    double nx = 0.0;
    double ny = 0.0;

    double u_hat = 0.0;
    double v_hat = 0.0;

    double uuh_hat = 0.0;
    double vvh_hat = 0.0;
    double uvh_hat = 0.0;
    double pe_hat  = 0.0;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_state.h_hat_at_gp[gp] = edge_state.ze_hat_at_gp[gp] + boundary.bath_at_gp[gp];

        nx = edge_bound.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.surface_normal[gp][GlobalCoord::y];

        u_hat = edge_state.qx_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];
        v_hat = edge_state.qy_hat_at_gp[gp] / edge_state.h_hat_at_gp[gp];

        uuh_hat = u_hat * edge_state.qx_hat_at_gp[gp];
        vvh_hat = v_hat * edge_state.qy_hat_at_gp[gp];
        uvh_hat = u_hat * edge_state.qy_hat_at_gp[gp];
        pe_hat  = Global::g * (0.5 * std::pow(edge_state.ze_hat_at_gp[gp], 2) +
                              edge_state.ze_hat_at_gp[gp] * boundary.bath_at_gp[gp]);

        // Add F_hat * n terms
        boundary.ze_numerical_flux_at_gp[gp] = edge_state.qx_hat_at_gp[gp] * nx + edge_state.qy_hat_at_gp[gp] * ny;
        boundary.qx_numerical_flux_at_gp[gp] = (uuh_hat + pe_hat) * nx + uvh_hat * ny;
        boundary.qy_numerical_flux_at_gp[gp] = uvh_hat * nx + (vvh_hat + pe_hat) * ny;
    }

    // Add tau * diff_q terms
    add_flux_tau_diff_q_LF(edge_bound.surface_normal,
                           boundary.ze_at_gp,
                           boundary.qx_at_gp,
                           boundary.qy_at_gp,
                           edge_state.ze_hat_at_gp,
                           edge_state.qx_hat_at_gp,
                           edge_state.qy_hat_at_gp,
                           edge_state.h_hat_at_gp,
                           boundary.ze_numerical_flux_at_gp,
                           boundary.qx_numerical_flux_at_gp,
                           boundary.qy_numerical_flux_at_gp);
}
}
}

#endif