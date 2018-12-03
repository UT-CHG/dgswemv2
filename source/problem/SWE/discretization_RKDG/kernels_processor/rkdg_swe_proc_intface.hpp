#ifndef RKDG_SWE_PROC_INTFACE_HPP
#define RKDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace RKDG {
template <typename InterfaceType>
void Problem::interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    auto& wd_state_in = intface.data_in.wet_dry_state;
    auto& wd_state_ex = intface.data_ex.wet_dry_state;

    if (wd_state_in.wet || wd_state_ex.wet) {
        const uint stage = stepper.GetStage();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        for ( uint var = 0; var < SWE::n_variables; ++var) {
            row(boundary_in.q_at_gp,var) = intface.ComputeUgpIN(state_in.q[var]);
            row(boundary_ex.q_at_gp,var) = intface.ComputeUgpEX(state_ex.q[var]);
        }
        intface.specialization.ComputeFlux(intface);

        // now compute contributions to the righthand side
        for ( uint var = 0; var < SWE::n_variables; ++var) {
            state_in.rhs[var] -= intface.IntegrationPhiIN(row(boundary_in.F_hat_at_gp,var));

            state_ex.rhs[var] -= intface.IntegrationPhiEX(row(boundary_ex.F_hat_at_gp,var));
        }
    }
}
}
}

#endif
