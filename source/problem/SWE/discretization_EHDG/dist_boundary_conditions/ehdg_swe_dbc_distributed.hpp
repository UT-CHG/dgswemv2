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
    auto& edge_state  = edge_dbound.edge_data.edge_state;
    auto& edge_global = edge_dbound.edge_data.edge_global;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    double ze_ex, qx_ex, qy_ex;
    double ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(
            gp, ze_ex, qx_ex, qy_ex, ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex);

        // Add F_in * n_in terms
        edge_global.ze_rhs_kernel_at_gp[gp] = -boundary.ze_flux_dot_n_at_gp[gp];
        edge_global.qx_rhs_kernel_at_gp[gp] = -boundary.qx_flux_dot_n_at_gp[gp];
        edge_global.qy_rhs_kernel_at_gp[gp] = -boundary.qy_flux_dot_n_at_gp[gp];

        // Add F_ex * n_ex terms
        edge_global.ze_rhs_kernel_at_gp[gp] += -ze_flux_dot_n_ex;
        edge_global.qx_rhs_kernel_at_gp[gp] += -qx_flux_dot_n_ex;
        edge_global.qy_rhs_kernel_at_gp[gp] += -qy_flux_dot_n_ex;
    }

    // Add tau terms
    add_kernel_tau_terms_dbound_LF(edge_dbound);
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        boundary.ze_numerical_flux_at_gp[gp] = boundary.ze_flux_dot_n_at_gp[gp];
        boundary.qx_numerical_flux_at_gp[gp] = boundary.qx_flux_dot_n_at_gp[gp];
        boundary.qy_numerical_flux_at_gp[gp] = boundary.qy_flux_dot_n_at_gp[gp];
    }

    // Add tau terms
    add_flux_tau_terms_bound_LF(edge_dbound);
}
}
}
}

#endif