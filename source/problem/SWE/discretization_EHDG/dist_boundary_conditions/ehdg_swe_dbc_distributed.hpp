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
}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_dbound.boundary.surface_normal;

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);
    SWE::get_dtau_dze_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dze);
    SWE::get_dtau_dqx_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqx);
    SWE::get_dtau_dqy_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqy);

    /* Get message from ex */

    std::vector<double> message;

    uint ngp = edge_dbound.boundary.data.get_ngp_boundary(edge_dbound.boundary.bound_id);

    message.resize(SWE::n_variables * ngp);

    edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);

    uint gp_ex;
    for (uint gp = 0; gp < ngp; ++gp) {
        gp_ex = ngp - gp - 1;

        for (uint var = 0; var < SWE::n_variables; ++var) {
            this->q_ex(var, gp_ex) = message[SWE::n_variables * gp + var];
        }
    }

    StatVector<double, SWE::n_variables> del_q;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

    for (uint gp = 0; gp < ngp; ++gp) {
        del_q = column(boundary.q_at_gp, gp) + column(this->q_ex, gp) - 2.0 * column(edge_internal.q_hat_at_gp, gp);

        column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
        column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
        column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

        column(edge_internal.rhs_global_kernel_at_gp, gp) = edge_internal.tau[gp] * del_q;

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) =
            flatten<double>(dtau_delq - 2.0 * edge_internal.tau[gp]);
    }

    // This whole calculation is pointless since in fact q_hat = 0.5 * (q_in + q_ex)
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    /* Compute trace flux */

    auto nx = row(edge_dbound.boundary.surface_normal, GlobalCoord::x);
    auto ny = row(edge_dbound.boundary.surface_normal, GlobalCoord::y);

    auto u = vec_cw_div(row(edge_internal.q_hat_at_gp, SWE::Variables::qx),
                        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(edge_internal.q_hat_at_gp, SWE::Variables::qy),
                        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(edge_internal.q_hat_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(edge_internal.q_hat_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(edge_internal.q_hat_at_gp, SWE::Variables::qy));
    auto pe  = Global::g * (0.5 * vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::ze),
                                             row(edge_internal.q_hat_at_gp, SWE::Variables::ze)) +
                           vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::ze),
                                       row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary.F_hat_at_gp, SWE::Variables::ze) =
        vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::qx), nx) +
        vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::qy), ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

    /* Add stabilization parameter terms */

    SWE::get_tau_LF(
        edge_internal.q_hat_at_gp, edge_internal.aux_hat_at_gp, edge_dbound.boundary.surface_normal, edge_internal.tau);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        column(boundary.F_hat_at_gp, gp) +=
            edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
    }
}
}
}
}

#endif