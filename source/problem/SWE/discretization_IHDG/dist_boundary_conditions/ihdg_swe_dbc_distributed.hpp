#ifndef IHDG_SWE_DBC_DISTRIBUTED_HPP
#define IHDG_SWE_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
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
    void Initialize(DistributedBoundaryType& dbound) {} /*nothing to initialize*/

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

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