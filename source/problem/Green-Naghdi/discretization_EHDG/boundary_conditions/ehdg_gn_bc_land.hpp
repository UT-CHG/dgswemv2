#ifndef EHDG_GN_BC_LAND_HPP
#define EHDG_GN_BC_LAND_HPP

namespace GN {
namespace EHDG {
namespace BC {
class Land : public SWE_SIM::BC::Land {
  public:
    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;
    auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    const auto nx = row(edge_bound.boundary.surface_normal, GlobalCoord::x);
    const auto ny = row(edge_bound.boundary.surface_normal, GlobalCoord::y);

    set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -1.0);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -1.0);

    set_constant(boundary.w1_hat_w1_kernel_at_gp, 0.0);
    set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), 1.0);
    set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), 1.0);
    row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx) -= vec_cw_mult(nx, nx);
    row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xy) -= vec_cw_mult(nx, ny);
    row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yx) -= vec_cw_mult(ny, nx);
    row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy) -= vec_cw_mult(ny, ny);
}
}
}
}

#endif