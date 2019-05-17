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

    // Our definition of numerical flux implies q_hat = 0.5 * (q_in + q_ex)
    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        row(edge_internal.q_hat_at_gp,var) = (boundary.q_at_gp[var] + row(q_ex, var)) / 2.0;
    }

    edge_state.q_hat = edge_dbound.L2Projection(edge_internal.q_hat_at_gp);
}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    boundary.delta_global_kernel_at_gp          = boundary.dF_hat_dq_at_gp;
    edge_internal.delta_hat_global_kernel_at_gp = boundary.dF_hat_dq_hat_at_gp;
    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        row(edge_internal.rhs_global_kernel_at_gp,var)       = boundary.F_hat_at_gp[var];
    }
}
}
}
}

#endif