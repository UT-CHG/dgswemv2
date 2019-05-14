#ifndef EHDG_SWE_PROC_INTFACE_HPP
#define EHDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace EHDG {
template <typename InterfaceType>
void Problem::global_interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    if (intface.data_in.wet_dry_state.wet || intface.data_ex.wet_dry_state.wet) {
        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            boundary_in.q_at_gp[var] = intface.ComputeUgpIN(state_in.q[var]);
            boundary_ex.q_at_gp[var] = intface.ComputeUgpEX(state_ex.q[var]);
        }
    }
}

template <typename InterfaceType>
void Problem::local_interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    if (intface.data_in.wet_dry_state.wet || intface.data_ex.wet_dry_state.wet) {
        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        // now compute contributions to the righthand side
        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            state_in.rhs[var] -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp[var]);

            state_ex.rhs[var] -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp[var]);
        }
    }
}
}
}

#endif
