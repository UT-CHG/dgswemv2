#ifndef RKDG_SWE_PROC_DBOUND_HPP
#define RKDG_SWE_PROC_DBOUND_HPP

namespace SWE {
namespace RKDG {
template <typename DistributedBoundaryType>
void Problem::distributed_boundary_send_kernel(const ProblemStepperType& stepper, DistributedBoundaryType& dbound) {
    auto& state      = dbound.data.state[stepper.GetStage()];
    auto& boundary   = dbound.data.boundary[dbound.bound_id];
    boundary.q_at_gp = dbound.ComputeUgp(state.q);
    row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
        row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
    std::vector<double> message;
    message.reserve(1 + (SWE::n_variables + 1) * ngp);
    message.push_back(dbound.data.wet_dry_state.wet);
    for (uint gp = 0; gp < ngp; ++gp) {
        message.push_back(boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));
        for (uint var = 0; var < SWE::n_variables; ++var) {
            message.push_back(boundary.q_at_gp(var, gp));
        }
    }
    dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::bound_state, message);
}

template <typename DistributedBoundaryType>
void Problem::distributed_boundary_kernel(const ProblemStepperType& stepper, DistributedBoundaryType& dbound) {
    std::vector<double> message(1);
    dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::bound_state, message);
    bool wet_ex = (bool)message[0];
    if (dbound.data.wet_dry_state.wet || wet_ex) {
        auto& state    = dbound.data.state[stepper.GetStage()];
        auto& boundary = dbound.data.boundary[dbound.bound_id];
        dbound.boundary_condition.ComputeFlux(dbound);
        state.rhs -= dbound.IntegrationPhi(boundary.F_hat_at_gp);
    }
}
}
}

#endif
