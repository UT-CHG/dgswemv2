#ifndef EHDG_GN_BC_FLOW_HPP
#define EHDG_GN_BC_FLOW_HPP

namespace GN {
namespace EHDG {
namespace BC {
class Flow : public SWE_SIM::BC::Flow {
  public:
    Flow() = default;
    Flow(const std::vector<SWE::FlowNode>& flow_input) : SWE_SIM::BC::Flow(flow_input) {}

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename StepperType, typename EdgeBoundaryType>
void Flow::ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;
    auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -1.0);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -1.0);

    set_constant(boundary.w1_hat_w1_kernel_at_gp, 0.0);
    set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), 1.0);
    set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), 1.0);
}
}
}
}

#endif