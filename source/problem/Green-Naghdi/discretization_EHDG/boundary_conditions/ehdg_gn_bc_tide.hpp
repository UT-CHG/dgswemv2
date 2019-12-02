#ifndef EHDG_GN_BC_TIDE_HPP
#define EHDG_GN_BC_TIDE_HPP

namespace GN {
namespace EHDG {
namespace BC {
class Tide : public SWE_SIM::BC::Tide {
  public:
    Tide() = default;
    Tide(const std::vector<SWE::TideNode>& tide_input) : SWE_SIM::BC::Tide(tide_input) {}

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename StepperType, typename EdgeBoundaryType>
void Tide::ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;
    auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double tau = -20;

    set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -tau);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -tau);

    set_constant(boundary.w1_hat_w1_kernel_at_gp, 0.0);
    set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), tau);
    set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), tau);

    boundary.w1_hat_w2_kernel_at_gp = edge_bound.boundary.surface_normal;
}
}
}
}

#endif
