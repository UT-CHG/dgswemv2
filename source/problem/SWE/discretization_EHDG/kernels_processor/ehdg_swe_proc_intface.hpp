#ifndef EHDG_SWE_PROC_INTFACE_HPP
#define EHDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace EHDG {
template <typename InterfaceType>
void Problem::global_interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    if (intface.data_in.wet_dry_state.wet || intface.data_ex.wet_dry_state.wet) {
        auto& state_in      = intface.data_in.state[stepper.GetStage()];
        auto& state_ex      = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];
        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);
        row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);
    }
}

template <typename InterfaceType>
void Problem::local_interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    if (intface.data_in.wet_dry_state.wet || intface.data_ex.wet_dry_state.wet) {
        auto& state_in    = intface.data_in.state[stepper.GetStage()];
        auto& state_ex    = intface.data_ex.state[stepper.GetStage()];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];
        state_in.rhs -= intface.IntegrationPhiIN(boundary_in.F_hat_at_gp);
        state_ex.rhs -= intface.IntegrationPhiEX(boundary_ex.F_hat_at_gp);
    }
}
}
}

#endif
