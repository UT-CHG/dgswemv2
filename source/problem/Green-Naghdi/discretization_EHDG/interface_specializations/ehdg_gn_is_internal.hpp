#ifndef EHDG_GN_IS_INTERNAL_HPP
#define EHDG_GN_IS_INTERNAL_HPP

namespace GN {
namespace EHDG {
namespace ISP {
class Internal : public SWE_SIM::ISP::Internal {
  public:
    template <typename EdgeInterfaceType>
    void ComputeGlobalKernelsDC(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernelsDC(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;
    auto& boundary_in   = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex   = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double tau = -20;

    if ((edge_int.interface.data_in.wet_dry_state.wet /*&& edge_int.interface.data_in.source.dispersive_correction*/) &&
        (edge_int.interface.data_ex.wet_dry_state.wet /*&& edge_int.interface.data_ex.source.dispersive_correction*/)) {
        set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -2.0 * tau);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -2.0 * tau);

        set_constant(boundary_in.w1_hat_w1_kernel_at_gp, 0.0);
        set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), tau);
        set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), tau);
        boundary_ex.w1_hat_w1_kernel_at_gp = boundary_in.w1_hat_w1_kernel_at_gp;

        boundary_in.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_in;
        boundary_ex.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_ex;
    } else if (edge_int.interface.data_in.wet_dry_state
                   .wet /*&& edge_int.interface.data_in.source.dispersive_correction*/) {
        /*set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -tau);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -tau);

        set_constant(boundary_in.w1_hat_w1_kernel_at_gp, 0.0);
        set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), tau);
        set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), tau);

        boundary_in.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_in;*/

        const auto nx = row(edge_int.interface.surface_normal_in, GlobalCoord::x);
        const auto ny = row(edge_int.interface.surface_normal_in, GlobalCoord::y);

        set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -1.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -1.0);

        set_constant(boundary_in.w1_hat_w1_kernel_at_gp, 0.0);
        set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), 1.0);
        set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), 1.0);
        row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx) -= vec_cw_mult(nx, nx);
        row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xy) -= vec_cw_mult(nx, ny);
        row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yx) -= vec_cw_mult(ny, nx);
        row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy) -= vec_cw_mult(ny, ny);

        set_constant(boundary_in.w1_hat_w2_kernel_at_gp, 0.0);
    } else if (edge_int.interface.data_ex.wet_dry_state
                   .wet /*&& edge_int.interface.data_ex.source.dispersive_correction*/) {
        /*set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -tau);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -tau);

        set_constant(boundary_ex.w1_hat_w1_kernel_at_gp, 0.0);
        set_constant(row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), tau);
        set_constant(row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), tau);

        boundary_ex.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_ex;*/

        const auto nx = row(edge_int.interface.surface_normal_ex, GlobalCoord::x);
        const auto ny = row(edge_int.interface.surface_normal_ex, GlobalCoord::y);

        set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -1.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -1.0);

        set_constant(boundary_ex.w1_hat_w1_kernel_at_gp, 0.0);
        set_constant(row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), 1.0);
        set_constant(row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), 1.0);
        row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx) -= vec_cw_mult(nx, nx);
        row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xy) -= vec_cw_mult(nx, ny);
        row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yx) -= vec_cw_mult(ny, nx);
        row(boundary_ex.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy) -= vec_cw_mult(ny, ny);

        set_constant(boundary_ex.w1_hat_w2_kernel_at_gp, 0.0);
    }
}
}
}
}

#endif