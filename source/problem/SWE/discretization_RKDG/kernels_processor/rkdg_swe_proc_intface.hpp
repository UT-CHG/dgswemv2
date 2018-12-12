#ifndef RKDG_SWE_PROC_INTFACE_HPP
#define RKDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace RKDG {
template <typename InterfaceType>
void Problem::interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    thread_local std::array<DynRowVector<double>, SWE::n_variables> q_in_at_gp;
    thread_local std::array<DynRowVector<double>, SWE::n_variables> q_ex_at_gp;

    auto& wd_state_in = intface.data_in.wet_dry_state;
    auto& wd_state_ex = intface.data_ex.wet_dry_state;

    if (wd_state_in.wet || wd_state_ex.wet) {
        const uint stage = stepper.GetStage();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        q_in_at_gp[SWE::Variables::ze] = intface.ComputeUgpIN(state_in.q[SWE::Variables::ze]);
        q_in_at_gp[SWE::Variables::qx] = intface.ComputeUgpIN(state_in.q[SWE::Variables::qx]);
        q_in_at_gp[SWE::Variables::qy] = intface.ComputeUgpIN(state_in.q[SWE::Variables::qy]);

        q_ex_at_gp[SWE::Variables::ze] = intface.ComputeUgpEX(state_ex.q[SWE::Variables::ze]);
        q_ex_at_gp[SWE::Variables::qx] = intface.ComputeUgpEX(state_ex.q[SWE::Variables::qx]);
        q_ex_at_gp[SWE::Variables::qy] = intface.ComputeUgpEX(state_ex.q[SWE::Variables::qy]);

        intface.specialization.ComputeFlux(intface,
                                           q_in_at_gp,
                                           q_ex_at_gp);

        // now compute contributions to the righthand side
        state_in.rhs[SWE::Variables::ze] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[SWE::Variables::ze]);
        state_in.rhs[SWE::Variables::qx] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[SWE::Variables::qx]);
        state_in.rhs[SWE::Variables::qy] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[SWE::Variables::qy]);

        state_ex.rhs[SWE::Variables::ze] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[SWE::Variables::ze]);
        state_ex.rhs[SWE::Variables::qx] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[SWE::Variables::qx]);
        state_ex.rhs[SWE::Variables::qy] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[SWE::Variables::qy]);

    }
}
}
}

#endif
