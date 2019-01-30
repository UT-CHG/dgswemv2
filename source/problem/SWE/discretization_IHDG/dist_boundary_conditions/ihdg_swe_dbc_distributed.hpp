#ifndef IHDG_SWE_DBC_DISTRIBUTED_HPP
#define IHDG_SWE_DBC_DISTRIBUTED_HPP

#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace IHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeInitTrace(EdgeDistributedType& edge_dbound, const HybMatrix<double, SWE::n_variables>& q_ex);

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {}

template <typename EdgeDistributedType>
void Distributed::ComputeInitTrace(EdgeDistributedType& edge_dbound, const HybMatrix<double, SWE::n_variables>& q_ex) {
    auto& dbound = edge_dbound.boundary;

    auto& boundary = dbound.data.boundary[dbound.bound_id];

    auto& edge_state    = edge_dbound.edge_data.edge_state;
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_dbound.boundary.surface_normal;

    set_constant(edge_state.q_hat, 0.0);

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        edge_internal.q_hat_at_gp = edge_dbound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);
        SWE::get_dtau_dze_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dze);
        SWE::get_dtau_dqx_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqx);
        SWE::get_dtau_dqy_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqy);

        StatVector<double, SWE::n_variables> del_q;
        StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

        for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
            del_q = column(boundary.q_at_gp, gp) + column(q_ex, gp) - 2.0 * column(q_hat_at_gp, gp);

            column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
            column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
            column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

            column(edge_internal.rhs_global_kernel_at_gp, gp) = edge_internal.tau[gp] * del_q;

            column(edge_internal.delta_hat_global_kernel_at_gp, gp) =
                flatten<double>(dtau_delq - 2.0 * edge_internal.tau[gp]);
        }

        for (uint dof_i = 0; dof_i < edge_dbound.edge_data.get_ndof(); ++dof_i) {
            for (uint dof_j = 0; dof_j < edge_dbound.edge_data.get_ndof(); ++dof_j) {
                submatrix(edge_internal.delta_hat_global,
                          SWE::n_variables * dof_i,
                          SWE::n_variables * dof_j,
                          SWE::n_variables,
                          SWE::n_variables) =
                    reshape<double, SWE::n_variables>(
                        edge_dbound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
            }

            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
                -edge_dbound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }

        solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

        edge_state.q_hat += reshape<double, SWE::n_variables, SO::ColumnMajor>(edge_internal.rhs_global,
                                                                               edge_dbound.edge_data.get_ndof());

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-12) {
            break;
        }
    }
}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    boundary.delta_global_kernel_at_gp          = boundary.dF_hat_dq_at_gp;
    edge_internal.delta_hat_global_kernel_at_gp = boundary.dF_hat_dq_hat_at_gp;
    edge_internal.rhs_global_kernel_at_gp       = boundary.F_hat_at_gp;
}
}
}
}

#endif