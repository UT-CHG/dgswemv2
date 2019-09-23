#ifndef EHDG_GN_PROC_EDGE_INTFACE_HPP
#define EHDG_GN_PROC_EDGE_INTFACE_HPP

namespace GN {
namespace EHDG {
template <typename EdgeInterfaceType>
void Problem::local_dc_edge_interface_kernel(const ESSPRKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;
    auto& internal_in = edge_int.interface.data_in.internal;
    auto& internal_ex = edge_int.interface.data_ex.internal;
    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double tau = -20;  // hardcode the tau value here

    // at this point h_at_gp
    // has been calculated in derivatives kernel
    // set h_hat as average of states
    const uint ngp = edge_int.edge_data.get_ngp();
    for (uint gp = 0; gp < ngp; ++gp) {
        const uint gp_ex = ngp - gp - 1;
        edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
            (boundary_in.aux_at_gp(SWE::Auxiliaries::h, gp) + boundary_ex.aux_at_gp(SWE::Auxiliaries::h, gp_ex)) / 2.0;
    }

    const auto h_hat = row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h);
    const auto bx = row(boundary_in.dbath_hat_at_gp, GlobalCoord::x);
    const auto by = row(boundary_in.dbath_hat_at_gp, GlobalCoord::y);
    const auto nx = row(edge_int.interface.surface_normal_in, GlobalCoord::x);
    const auto ny = row(edge_int.interface.surface_normal_in, GlobalCoord::y);

    set_constant(boundary_in.w1_w1_kernel_at_gp, 0.0);
    set_constant(row(boundary_in.w1_w1_kernel_at_gp, RowMajTrans2D::xx), -NDParameters::alpha / 3.0 * tau);
    set_constant(row(boundary_in.w1_w1_kernel_at_gp, RowMajTrans2D::yy), -NDParameters::alpha / 3.0 * tau);
    boundary_ex.w1_w1_kernel_at_gp = boundary_in.w1_w1_kernel_at_gp;

    set_constant(boundary_in.w1_w1_hat_kernel_at_gp, 0.0);
    set_constant(row(boundary_in.w1_w1_hat_kernel_at_gp, RowMajTrans2D::xx), NDParameters::alpha / 3.0 * tau);
    set_constant(row(boundary_in.w1_w1_hat_kernel_at_gp, RowMajTrans2D::yy), NDParameters::alpha / 3.0 * tau);
    boundary_ex.w1_w1_hat_kernel_at_gp = boundary_in.w1_w1_hat_kernel_at_gp;

    row(boundary_in.w1_w1_hat_kernel_at_gp, RowMajTrans2D::xx) +=
            NDParameters::alpha / 2.0 * vec_cw_mult(vec_cw_mult(bx, nx), h_hat);
    row(boundary_in.w1_w1_hat_kernel_at_gp, RowMajTrans2D::xy) +=
            NDParameters::alpha / 2.0 * vec_cw_mult(vec_cw_mult(by, nx), h_hat);
    row(boundary_in.w1_w1_hat_kernel_at_gp, RowMajTrans2D::yx) +=
            NDParameters::alpha / 2.0 * vec_cw_mult(vec_cw_mult(bx, ny), h_hat);
    row(boundary_in.w1_w1_hat_kernel_at_gp, RowMajTrans2D::yy) +=
            NDParameters::alpha / 2.0 * vec_cw_mult(vec_cw_mult(by, ny), h_hat);

    row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::x) = -vec_cw_div(nx, h_hat);
    row(boundary_in.w2_w1_hat_kernel_at_gp, GlobalCoord::y) = -vec_cw_div(ny, h_hat);

    for (uint gp = 0; gp < ngp; ++gp) {
        const uint gp_ex = ngp - gp - 1;
        const double h_hat = edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);
        const double bx = boundary_ex.dbath_hat_at_gp(GlobalCoord::x, gp_ex);
        const double by = boundary_ex.dbath_hat_at_gp(GlobalCoord::y, gp_ex);
        const double nx = edge_int.interface.surface_normal_ex(GlobalCoord::x, gp_ex);
        const double ny = edge_int.interface.surface_normal_ex(GlobalCoord::y, gp_ex);

        boundary_ex.w1_w1_hat_kernel_at_gp(RowMajTrans2D::xx, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat * bx * nx;
        boundary_ex.w1_w1_hat_kernel_at_gp(RowMajTrans2D::xy, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat * by * nx;
        boundary_ex.w1_w1_hat_kernel_at_gp(RowMajTrans2D::yx, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat * bx * ny;
        boundary_ex.w1_w1_hat_kernel_at_gp(RowMajTrans2D::yy, gp_ex) +=
            NDParameters::alpha / 2.0 * h_hat * by * ny;

        boundary_ex.w2_w1_hat_kernel_at_gp(GlobalCoord::x, gp_ex) = -nx / h_hat;
        boundary_ex.w2_w1_hat_kernel_at_gp(GlobalCoord::y, gp_ex) = -ny / h_hat;
    }

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
void Problem::global_dc_edge_interface_kernel(const ESSPRKStepper& stepper, EdgeInterfaceType& edge_int) {
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