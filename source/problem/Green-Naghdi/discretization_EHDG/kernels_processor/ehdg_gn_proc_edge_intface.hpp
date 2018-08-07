#ifndef EHDG_GN_PROC_EDGE_INTFACE_HPP
#define EHDG_GN_PROC_EDGE_INTFACE_HPP

namespace GN {
namespace EHDG {
template <typename EdgeInterfaceType>
void Problem::global_swe_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    /* Take average of in/ex state as initial trace state */
    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(edge_internal.q_init_at_gp, gp) =
            (column(boundary_in.q_at_gp, gp) + column(boundary_ex.q_at_gp, gp_ex)) / 2.0;
    }

    edge_state.q_hat = edge_int.L2Projection(edge_internal.q_init_at_gp);

    /* Newton-Raphson iterator */

    uint iter = 0;
    while (true) {
        iter++;

        Problem::global_swe_edge_interface_iteration(stepper, edge_int);

        if (iter == 100) {
            break;
        }

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-8) {
            break;
        }
    }

    /* Compute Numerical Flux */

    edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, GN::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, GN::Variables::ze) + row(boundary_in.aux_at_gp, GN::Auxiliaries::bath);

    edge_int.interface.specialization.ComputeNumericalFluxSWE(edge_int);
}

template <typename EdgeInterfaceType>
void Problem::global_swe_edge_interface_iteration(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];

    edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, GN::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, GN::Variables::ze) + row(boundary_in.aux_at_gp, GN::Auxiliaries::bath);

    /* Assemble global kernels */

    edge_int.interface.specialization.ComputeGlobalKernelsSWE(edge_int);

    /* Assemble global system */

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      GN::n_variables * dof_i,
                      GN::n_variables * dof_j,
                      GN::n_variables,
                      GN::n_variables) =
                reshape<double, GN::n_variables>(
                    edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, GN::n_variables * dof_i, GN::n_variables) =
            -edge_int.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    /* Solve global system for delta_q */

    solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

    /* Increment q */

    for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
        edge_state.q_hat(GN::Variables::ze, dof) += edge_internal.rhs_global[3 * dof];
        edge_state.q_hat(GN::Variables::qx, dof) += edge_internal.rhs_global[3 * dof + 1];
        edge_state.q_hat(GN::Variables::qy, dof) += edge_internal.rhs_global[3 * dof + 2];
    }
}

template <typename EdgeInterfaceType>
void Problem::local_dc_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& internal_in = edge_int.interface.data_in.internal;
    auto& internal_ex = edge_int.interface.data_ex.internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    // at this point h_hat_at_gp
    // has been computed in swe edge_boundary kernel

    double tau = -20.0;  // hardcode the tau value here

    // set kernels up
    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(boundary_in.w1_w1_kernel_at_gp, gp) =
            -NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);
        column(boundary_in.w1_w1_hat_kernel_at_gp, gp) =
            NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);

        column(boundary_ex.w1_w1_kernel_at_gp, gp_ex) =
            -NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);
        column(boundary_ex.w1_w1_hat_kernel_at_gp, gp_ex) =
            NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);
    }

    row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::x) =
        -cwise_division(row(edge_int.interface.surface_normal_in, GlobalCoord::x),
                        row(edge_internal.aux_hat_at_gp, GN::Auxiliaries::h));
    row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::y) =
        -cwise_division(row(edge_int.interface.surface_normal_in, GlobalCoord::y),
                        row(edge_internal.aux_hat_at_gp, GN::Auxiliaries::h));

    row(boundary_ex.w2_w1_hat_kernel_at_gp, GlobalCoord::x) =
        -cwise_division(row(edge_int.interface.surface_normal_ex, GlobalCoord::x),
                        row(edge_internal.aux_hat_at_gp, GN::Auxiliaries::h));
    row(boundary_ex.w2_w1_hat_kernel_at_gp, GlobalCoord::y) =
        -cwise_division(row(edge_int.interface.surface_normal_ex, GlobalCoord::y),
                        row(edge_internal.aux_hat_at_gp, GN::Auxiliaries::h));

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); dof_j++) {
            submatrix(internal_in.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) +=
                reshape<double, GN::n_dimensions>(
                    edge_int.interface.IntegrationPhiPhiIN(dof_i, dof_j, boundary_in.w1_w1_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); dof_j++) {
            submatrix(boundary_in.w1_w1_hat,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_int.IntegrationPhiLambdaIN(dof_i, dof_j, boundary_in.w1_w1_hat_kernel_at_gp));

            boundary_in.w2_w1_hat(dof_i, GN::n_dimensions * dof_j + GlobalCoord::x) =
                edge_int.IntegrationPhiLambdaIN(dof_i, dof_j, row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::x));

            boundary_in.w2_w1_hat(dof_i, GN::n_dimensions * dof_j + GlobalCoord::y) =
                edge_int.IntegrationPhiLambdaIN(dof_i, dof_j, row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::y));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); dof_j++) {
            submatrix(internal_ex.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) +=
                reshape<double, GN::n_dimensions>(
                    edge_int.interface.IntegrationPhiPhiEX(dof_i, dof_j, boundary_ex.w1_w1_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); dof_j++) {
            submatrix(boundary_ex.w1_w1_hat,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_int.IntegrationPhiLambdaEX(dof_i, dof_j, boundary_ex.w1_w1_hat_kernel_at_gp));

            boundary_ex.w2_w1_hat(dof_i, GN::n_dimensions * dof_j + GlobalCoord::x) =
                edge_int.IntegrationPhiLambdaEX(dof_i, dof_j, row(boundary_ex.w2_w1_hat_kernel_at_gp, GlobalCoord::x));

            boundary_ex.w2_w1_hat(dof_i, GN::n_dimensions * dof_j + GlobalCoord::y) =
                edge_int.IntegrationPhiLambdaEX(dof_i, dof_j, row(boundary_ex.w2_w1_hat_kernel_at_gp, GlobalCoord::y));
        }
    }
}

template <typename EdgeInterfaceType>
void Problem::global_dc_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.interface.specialization.ComputeGlobalKernelsDC(edge_int);

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); dof_j++) {
            submatrix(edge_internal.w1_hat_w1_hat,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.w1_hat_w1_hat_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); dof_j++) {
            submatrix(boundary_in.w1_hat_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_int.IntegrationPhiLambdaIN(dof_j, dof_i, boundary_in.w1_hat_w1_kernel_at_gp));

            boundary_in.w1_hat_w2(GN::n_dimensions * dof_i + GlobalCoord::x, dof_j) =
                edge_int.IntegrationPhiLambdaIN(dof_j, dof_i, row(boundary_in.w1_hat_w2_kernel_at_gp, GlobalCoord::x));

            boundary_in.w1_hat_w2(GN::n_dimensions * dof_i + GlobalCoord::y, dof_j) =
                edge_int.IntegrationPhiLambdaIN(dof_j, dof_i, row(boundary_in.w1_hat_w2_kernel_at_gp, GlobalCoord::y));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); dof_j++) {
            submatrix(boundary_ex.w1_hat_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_int.IntegrationPhiLambdaEX(dof_j, dof_i, boundary_ex.w1_hat_w1_kernel_at_gp));

            boundary_ex.w1_hat_w2(GN::n_dimensions * dof_i + GlobalCoord::x, dof_j) =
                edge_int.IntegrationPhiLambdaEX(dof_j, dof_i, row(boundary_ex.w1_hat_w2_kernel_at_gp, GlobalCoord::x));

            boundary_ex.w1_hat_w2(GN::n_dimensions * dof_i + GlobalCoord::y, dof_j) =
                edge_int.IntegrationPhiLambdaEX(dof_j, dof_i, row(boundary_ex.w1_hat_w2_kernel_at_gp, GlobalCoord::y));
        }
    }
}
}
}

#endif