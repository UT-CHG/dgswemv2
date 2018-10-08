#ifndef EHDG_GN_PROC_EDGE_INTFACE_HPP
#define EHDG_GN_PROC_EDGE_INTFACE_HPP

namespace GN {
namespace EHDG {
template <typename EdgeInterfaceType>
void Problem::local_dc_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& internal_in = edge_int.interface.data_in.internal;
    auto& internal_ex = edge_int.interface.data_ex.internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    // at this point h_at_gp
    // has been calculated in derivatives kernel

    // set h_hat as average of states
    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
            (boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp) + boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex)) / 2.0;
    }

    double tau = -20;  // hardcode the tau value here

    // set kernels up
    double h_hat_in, h_hat_ex;
    double bx_in, by_in, bx_ex, by_ex;
    double nx_in, ny_in, nx_ex, ny_ex;

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        h_hat_in = edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        h_hat_ex = edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp_ex);

        bx_in = boundary_in.dbath_hat_at_gp(GlobalCoord::x, gp);
        by_in = boundary_in.dbath_hat_at_gp(GlobalCoord::y, gp);
        nx_in = edge_int.interface.surface_normal_in(GlobalCoord::x, gp);
        ny_in = edge_int.interface.surface_normal_in(GlobalCoord::y, gp);

        bx_ex = boundary_ex.dbath_hat_at_gp(GlobalCoord::x, gp_ex);
        by_ex = boundary_ex.dbath_hat_at_gp(GlobalCoord::y, gp_ex);
        nx_ex = edge_int.interface.surface_normal_ex(GlobalCoord::x, gp_ex);
        ny_ex = edge_int.interface.surface_normal_ex(GlobalCoord::y, gp_ex);

        column(boundary_in.w1_w1_kernel_at_gp, gp) =
            -NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);
        column(boundary_in.w1_w1_hat_kernel_at_gp, gp) =
            NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);

        boundary_in.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::x, gp) +=
            NDParameters::alpha / 2.0 * h_hat_in * bx_in * nx_in;
        boundary_in.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::y, gp) +=
            NDParameters::alpha / 2.0 * h_hat_in * by_in * nx_in;
        boundary_in.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::x, gp) +=
            NDParameters::alpha / 2.0 * h_hat_in * bx_in * ny_in;
        boundary_in.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::y, gp) +=
            NDParameters::alpha / 2.0 * h_hat_in * by_in * ny_in;

        column(boundary_ex.w1_w1_kernel_at_gp, gp_ex) =
            -NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);
        column(boundary_ex.w1_w1_hat_kernel_at_gp, gp_ex) =
            NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);

        boundary_ex.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::x, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat_ex * bx_ex * nx_ex;
        boundary_ex.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::y, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat_ex * by_ex * nx_ex;
        boundary_ex.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::x, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat_ex * bx_ex * ny_ex;
        boundary_ex.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::y, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat_ex * by_ex * ny_ex;
    }

    row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::x) =
        -vec_cw_div(row(edge_int.interface.surface_normal_in, GlobalCoord::x),
                    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));
    row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::y) =
        -vec_cw_div(row(edge_int.interface.surface_normal_in, GlobalCoord::y),
                    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));

    row(boundary_ex.w2_w1_hat_kernel_at_gp, GlobalCoord::x) =
        -vec_cw_div(row(edge_int.interface.surface_normal_ex, GlobalCoord::x),
                    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));
    row(boundary_ex.w2_w1_hat_kernel_at_gp, GlobalCoord::y) =
        -vec_cw_div(row(edge_int.interface.surface_normal_ex, GlobalCoord::y),
                    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); ++dof_j) {
            submatrix(internal_in.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) +=
                reshape<double, GN::n_dimensions>(
                    edge_int.interface.IntegrationPhiPhiIN(dof_i, dof_j, boundary_in.w1_w1_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
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

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); ++dof_j) {
            submatrix(internal_ex.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) +=
                reshape<double, GN::n_dimensions>(
                    edge_int.interface.IntegrationPhiPhiEX(dof_i, dof_j, boundary_ex.w1_w1_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
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

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.w1_hat_w1_hat,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.w1_hat_w1_hat_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); ++dof_j) {
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

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); ++dof_j) {
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