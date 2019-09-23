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
    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double tau = -20;

    set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -2.0 * tau);
    set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -2.0 * tau);

    set_constant(boundary_in.w1_hat_w1_kernel_at_gp, 0.0);
    set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), tau);
    set_constant(row(boundary_in.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), tau);
    boundary_ex.w1_hat_w1_kernel_at_gp = boundary_in.w1_hat_w1_kernel_at_gp;

    boundary_in.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_in;
    boundary_ex.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_ex;
}
}
}
}

#endif