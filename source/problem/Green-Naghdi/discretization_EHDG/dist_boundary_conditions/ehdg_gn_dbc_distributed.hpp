#ifndef EHDG_GN_DBC_DISTRIBUTED_HPP
#define EHDG_GN_DBC_DISTRIBUTED_HPP

#include "communication/db_data_exchanger.hpp"

namespace GN {
namespace EHDG {
namespace DBC {
class Distributed : public SWE_SIM::DBC::Distributed {
  public:
    bool wet_neighbor = false;

    Distributed(const DBDataExchanger& exchanger) : SWE_SIM::DBC::Distributed(exchanger) {}

    template <typename EdgeDistributedType>
    void ComputeGlobalKernelsDC(EdgeDistributedType& edge_dbound);
};

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernelsDC(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;
    auto& boundary      = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    double tau = -20;

    if (wet_neighbor) {
        set_constant(edge_internal.w1_hat_w1_hat_kernel_at_gp, 0.0);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::xx), -tau);
        set_constant(row(edge_internal.w1_hat_w1_hat_kernel_at_gp, RowMajTrans2D::yy), -tau);

        set_constant(boundary.w1_hat_w1_kernel_at_gp, 0.0);
        set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::xx), tau);
        set_constant(row(boundary.w1_hat_w1_kernel_at_gp, RowMajTrans2D::yy), tau);

        boundary.w1_hat_w2_kernel_at_gp = edge_dbound.boundary.surface_normal;
    } else {
        const auto nx = row(edge_dbound.boundary.surface_normal, GlobalCoord::x);
        const auto ny = row(edge_dbound.boundary.surface_normal, GlobalCoord::y);

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

        set_constant(boundary.w1_hat_w2_kernel_at_gp, 0.0);
    }
}
}
}
}

#endif