#ifndef EHDG_GN_PROC_EDGE_BOUND_HPP
#define EHDG_GN_PROC_EDGE_BOUND_HPP

namespace GN {
namespace EHDG {
template <typename EdgeBoundaryType>
void Problem::local_dc_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& internal = edge_bound.boundary.data.internal;
    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    // at this point h_at_gp
    // has been calculated in derivatives kernel

    // set h_hat as internal state
    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) = row(boundary.aux_at_gp, SWE::Auxiliaries::h);

    double tau = -20;  // hardcode the tau value here

    // set kernels up
    double h_hat;
    double bx, by;
    double nx, ny;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        h_hat = edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp);

        bx = boundary.dbath_hat_at_gp(GlobalCoord::x, gp);
        by = boundary.dbath_hat_at_gp(GlobalCoord::y, gp);

        nx = edge_bound.boundary.surface_normal(GlobalCoord::x, gp);
        ny = edge_bound.boundary.surface_normal(GlobalCoord::y, gp);

        column(boundary.w1_w1_kernel_at_gp, gp) =
            -NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);

        column(boundary.w1_w1_hat_kernel_at_gp, gp) =
            NDParameters::alpha / 3.0 * tau * IdentityVector<double>(GN::n_dimensions);

        boundary.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::x, gp) +=
            NDParameters::alpha / 2.0 * h_hat * bx * nx;
        boundary.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::y, gp) +=
            NDParameters::alpha / 2.0 * h_hat * by * nx;
        boundary.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::x, gp) +=
            NDParameters::alpha / 2.0 * h_hat * bx * ny;
        boundary.w1_w1_hat_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::y, gp) +=
            NDParameters::alpha / 2.0 * h_hat * by * ny;
    }

    row(boundary.w2_w1_hat_kernel_at_gp, GlobalCoord::x) = -vec_cw_div(
        row(edge_bound.boundary.surface_normal, GlobalCoord::x), row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));
    row(boundary.w2_w1_hat_kernel_at_gp, GlobalCoord::y) = -vec_cw_div(
        row(edge_bound.boundary.surface_normal, GlobalCoord::y), row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(internal.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) +=
                reshape<double, GN::n_dimensions>(
                    edge_bound.boundary.IntegrationPhiPhi(dof_i, dof_j, boundary.w1_w1_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
            submatrix(boundary.w1_w1_hat,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_bound.IntegrationPhiLambda(dof_i, dof_j, boundary.w1_w1_hat_kernel_at_gp));

            boundary.w2_w1_hat(dof_i, GN::n_dimensions * dof_j + GlobalCoord::x) =
                edge_bound.IntegrationPhiLambda(dof_i, dof_j, row(boundary.w2_w1_hat_kernel_at_gp, GlobalCoord::x));

            boundary.w2_w1_hat(dof_i, GN::n_dimensions * dof_j + GlobalCoord::y) =
                edge_bound.IntegrationPhiLambda(dof_i, dof_j, row(boundary.w2_w1_hat_kernel_at_gp, GlobalCoord::y));
        }
    }
}

template <typename EdgeBoundaryType>
void Problem::global_dc_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_bound.boundary.boundary_condition.ComputeGlobalKernelsDC(stepper, edge_bound);

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.w1_hat_w1_hat,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.w1_hat_w1_hat_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(boundary.w1_hat_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(
                    edge_bound.IntegrationPhiLambda(dof_j, dof_i, boundary.w1_hat_w1_kernel_at_gp));
        }
    }
}
}
}

#endif