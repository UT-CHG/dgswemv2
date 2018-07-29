#ifndef EHDG_GN_DBC_DISTRIBUTED_HPP
#define EHDG_GN_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "ehdg_gn_dbc_db_data_exch.hpp"

namespace GN {
namespace EHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

    HybMatrix<double, GN::n_variables> q_ex;
    HybMatrix<double, GN::n_variables> Fn_ex;

  public:
    Distributed() = default;
    Distributed(const DBDataExchanger& exchanger);

    template <typename DistributedBoundaryType>
    void Initialize(DistributedBoundaryType& dbound);

    template <typename EdgeDistributedType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound);

    template <typename EdgeDistributedType>
    void ComputeNumericalFlux(EdgeDistributedType& edge_dbound);
};

Distributed::Distributed(const DBDataExchanger& exchanger) : exchanger(exchanger) {}

template <typename DistributedBoundaryType>
void Distributed::Initialize(DistributedBoundaryType& dbound) {
    uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    this->q_ex.resize(GN::n_variables, ngp);
    this->Fn_ex.resize(GN::n_variables, ngp);
}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    StatVector<double, GN::n_variables> q_ex;
    StatVector<double, GN::n_variables> Fn_ex;

    edge_internal.rhs_global_kernel_at_gp = boundary.Fn_at_gp + this->Fn_ex;

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