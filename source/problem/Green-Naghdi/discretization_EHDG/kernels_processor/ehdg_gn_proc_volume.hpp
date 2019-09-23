#ifndef EHDG_GN_PROC_VOLUME_HPP
#define EHDG_GN_PROC_VOLUME_HPP

namespace GN {
namespace EHDG {
template <typename ElementType>
void Problem::local_dc_volume_kernel(const ESSPRKStepper& stepper, ElementType& elt) {
    auto& internal = elt.data.internal;

    // at this point h_at_gp
    // has been calculated in derivatives kernel

    const auto h  = row(internal.aux_at_gp, SWE::Auxiliaries::h);
    const auto bx = row(internal.dbath_at_gp, GlobalCoord::x);
    const auto by = row(internal.dbath_at_gp, GlobalCoord::y);

    set_constant(internal.w1_w1_kernel_at_gp, 0.0);
    set_constant(row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xx), 1.0);
    set_constant(row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yy), 1.0);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xx) += NDParameters::alpha * vec_cw_mult(bx, bx);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xy) += NDParameters::alpha * vec_cw_mult(bx, by);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yx) += NDParameters::alpha * vec_cw_mult(by, bx);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yy) += NDParameters::alpha * vec_cw_mult(by, by);

    set_constant(internal.w1_w2_kernel_at_gp, -NDParameters::alpha / 3.0);
    row(internal.w2_w1_kernel_at_gp, GlobalCoord::x) = power(h, -1.0);
    row(internal.w2_w1_kernel_at_gp, GlobalCoord::y) = power(h, -1.0);
    internal.w2_w2_kernel_at_gp                      = power(h, -3.0);

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::x, dof_j) =
                elt.IntegrationPhiDPhi(dof_i, GlobalCoord::x, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::x));
            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::y, dof_j) =
                elt.IntegrationPhiDPhi(dof_i, GlobalCoord::y, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::y));
        }
    }

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            submatrix(internal.w1_w1,
                      GN::n_dimensions * dof_i,
                      GN::n_dimensions * dof_j,
                      GN::n_dimensions,
                      GN::n_dimensions) =
                reshape<double, GN::n_dimensions>(elt.IntegrationPhiPhi(dof_i, dof_j, internal.w1_w1_kernel_at_gp));
            internal.w2_w1(dof_i, GN::n_dimensions * dof_j + GlobalCoord::x) =
                elt.IntegrationPhiDPhi(dof_j, GlobalCoord::x, dof_i, row(internal.w2_w1_kernel_at_gp, GlobalCoord::x));
            internal.w2_w1(dof_i, GN::n_dimensions * dof_j + GlobalCoord::y) =
                elt.IntegrationPhiDPhi(dof_j, GlobalCoord::y, dof_i, row(internal.w2_w1_kernel_at_gp, GlobalCoord::y));
            internal.w2_w2(dof_i, dof_j) = elt.IntegrationPhiPhi(dof_i, dof_j, internal.w2_w2_kernel_at_gp);
        }
    }

    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xx) = -NDParameters::alpha / 2.0 * vec_cw_mult(h, bx);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xy) = -NDParameters::alpha / 2.0 * vec_cw_mult(h, by);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yx) = -NDParameters::alpha / 2.0 * vec_cw_mult(h, bx);
    row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yy) = -NDParameters::alpha / 2.0 * vec_cw_mult(h, by);

    row(internal.w1_w2_kernel_at_gp, GlobalCoord::x) = -NDParameters::alpha / 2.0 * vec_cw_div(bx, h);
    row(internal.w1_w2_kernel_at_gp, GlobalCoord::y) = -NDParameters::alpha / 2.0 * vec_cw_div(by, h);

    for (uint dof_i = 0; dof_i < elt.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < elt.data.get_ndof(); ++dof_j) {
            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::x, GN::n_dimensions * dof_j + GlobalCoord::x) +=
                elt.IntegrationPhiDPhi(
                    dof_j, GlobalCoord::x, dof_i, row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xx));
            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::x, GN::n_dimensions * dof_j + GlobalCoord::y) +=
                elt.IntegrationPhiDPhi(
                    dof_j, GlobalCoord::x, dof_i, row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::xy));
            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::y, GN::n_dimensions * dof_j + GlobalCoord::x) +=
                elt.IntegrationPhiDPhi(
                    dof_j, GlobalCoord::y, dof_i, row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yx));
            internal.w1_w1(GN::n_dimensions * dof_i + GlobalCoord::y, GN::n_dimensions * dof_j + GlobalCoord::y) +=
                elt.IntegrationPhiDPhi(
                    dof_j, GlobalCoord::y, dof_i, row(internal.w1_w1_kernel_at_gp, RowMajTrans2D::yy));
            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::x, dof_j) +=
                elt.IntegrationPhiPhi(dof_i, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::x));
            internal.w1_w2(GN::n_dimensions * dof_i + GlobalCoord::y, dof_j) +=
                elt.IntegrationPhiPhi(dof_i, dof_j, row(internal.w1_w2_kernel_at_gp, GlobalCoord::y));
        }
    }
}

template <typename ElementType>
void Problem::dispersive_correction_kernel(const ESSPRKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();
    auto& state      = elt.data.state[stage];
    auto& internal   = elt.data.internal;

    const auto h      = row(internal.aux_at_gp, SWE::Auxiliaries::h);
    const auto dze_dx = elt.ComputeUgp(row(state.dze, GlobalCoord::x));
    const auto dze_dy = elt.ComputeUgp(row(state.dze, GlobalCoord::y));

    row(internal.source_at_gp, SWE::Variables::qx) = Global::g / NDParameters::alpha * vec_cw_mult(dze_dx, h);
    row(internal.source_at_gp, SWE::Variables::qy) = Global::g / NDParameters::alpha * vec_cw_mult(dze_dy, h);
    row(internal.source_at_gp, SWE::Variables::qx) -= elt.ComputeUgp(row(state.w1, GlobalCoord::x));
    row(internal.source_at_gp, SWE::Variables::qy) -= elt.ComputeUgp(row(state.w1, GlobalCoord::y));

    set_constant(row(state.rhs, SWE::Variables::ze), 0);
    row(state.rhs, SWE::Variables::qx) = elt.IntegrationPhi(row(internal.source_at_gp, SWE::Variables::qx));
    row(state.rhs, SWE::Variables::qy) = elt.IntegrationPhi(row(internal.source_at_gp, SWE::Variables::qy));
}
}
}

#endif
