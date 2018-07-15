#ifndef EHDG_SWE_DBC_DISTRIBUTED_HPP
#define EHDG_SWE_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "ehdg_swe_dbc_db_data_exch.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound) {} /*nothing to initialize*/

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    StatVector<double, SWE::n_variables> q_ex;
    StatVector<double, SWE::n_variables> Fn_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(gp, q_ex, Fn_ex);

        edge_internal.rhs_global_kernel_at_gp[gp] = boundary.Fn_at_gp[gp];
        edge_internal.rhs_global_kernel_at_gp[gp] += Fn_ex;
    }

    // Add tau terms
    add_kernel_tau_terms_dbound_LF(edge_dbound);
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    boundary.F_hat_at_gp = boundary.Fn_at_gp;

    // Add tau terms
    add_F_hat_tau_terms_bound_LF(edge_dbound);
}
}
}
}

#endif