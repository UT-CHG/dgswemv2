#ifndef EHDG_SWE_DBC_DISTRIBUTED_HPP
#define EHDG_SWE_DBC_DISTRIBUTED_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "communication/db_data_exchanger.hpp"

namespace SWE {
namespace EHDG {
namespace DBC {
class Distributed {
  public:
    DBDataExchanger exchanger;

    HybMatrix<double, SWE::n_variables> q_ex;
    HybMatrix<double, SWE::n_variables> Fn_ex;

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
    this->q_ex.resize(SWE::n_variables, ngp);
    this->Fn_ex.resize(SWE::n_variables, ngp);
}

template <typename EdgeDistributedType>
void Distributed::ComputeGlobalKernels(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

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

    edge_internal.rhs_global_kernel_at_gp = boundary.Fn_at_gp + this->Fn_ex;

    // Add tau terms
    SWE::EHDG::add_kernel_tau_terms_dbound_LF(edge_dbound);
}

template <typename EdgeDistributedType>
void Distributed::ComputeNumericalFlux(EdgeDistributedType& edge_dbound) {
    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    boundary.F_hat_at_gp = boundary.Fn_at_gp;

    // Add tau terms
    SWE::EHDG::add_F_hat_tau_terms_bound_LF(edge_dbound);
}
}
}
}

#endif