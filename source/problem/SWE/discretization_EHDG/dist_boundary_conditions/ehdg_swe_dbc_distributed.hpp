#ifndef EHDG_SWE_DBC_DISTRIBUTED_HPP
#define EHDG_SWE_DBC_DISTRIBUTED_HPP

#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

  private:
    HybMatrix<double, SWE::n_variables> q_ex;
    HybMatrix<double, SWE::n_variables> Fn_ex;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(EdgeDistributedType& edge_dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {
    uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

    this->q_ex.resize(SWE::n_variables, ngp);
    this->Fn_ex.resize(SWE::n_variables, ngp);
}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    SWE::get_tau_LF(
        edge_internal.q_hat_at_gp, edge_internal.aux_hat_at_gp, edge_dbound.boundary.surface_normal, edge_internal.tau);
    SWE::get_dtau_dze_LF(edge_internal.q_hat_at_gp,
                         edge_internal.aux_hat_at_gp,
                         edge_dbound.boundary.surface_normal,
                         edge_internal.dtau_dze);
    SWE::get_dtau_dqx_LF(edge_internal.q_hat_at_gp,
                         edge_internal.aux_hat_at_gp,
                         edge_dbound.boundary.surface_normal,
                         edge_internal.dtau_dqx);
    SWE::get_dtau_dqy_LF(edge_internal.q_hat_at_gp,
                         edge_internal.aux_hat_at_gp,
                         edge_dbound.boundary.surface_normal,
                         edge_internal.dtau_dqy);

    std::vector<double> message;

    uint ngp = edge_dbound.boundary.data.get_ngp_boundary(edge_dbound.boundary.bound_id);

    message.resize(2 * SWE::n_variables * ngp);

    edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);

    uint gp_ex;
    for (uint gp = 0; gp < ngp; ++gp) {
        gp_ex = ngp - gp - 1;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            this->q_ex(var, gp_ex)  = message[SWE::n_variables * gp + var];
            this->Fn_ex(var, gp_ex) = message[SWE::n_variables * ngp + SWE::n_variables * gp + var];
        }
    }

    StatVector<double, SWE::n_variables> del_q;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

    for (uint gp = 0; gp < ngp; ++gp) {
        del_q = column(boundary.q_at_gp, gp) + column(this->q_ex, gp) - 2.0 * column(edge_internal.q_hat_at_gp, gp);

        column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
        column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
        column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            column(boundary.Fn_at_gp, gp) + column(this->Fn_ex, gp) + edge_internal.tau[gp] * del_q;

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) =
            flatten<double>(dtau_delq - 2.0 * edge_internal.tau[gp]);
    }
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    SWE::get_tau_LF(
        edge_internal.q_hat_at_gp, edge_internal.aux_hat_at_gp, edge_dbound.boundary.surface_normal, edge_internal.tau);

    boundary.F_hat_at_gp = boundary.Fn_at_gp;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        column(boundary.F_hat_at_gp, gp) +=
            edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
    }
}
}
}
}

#endif