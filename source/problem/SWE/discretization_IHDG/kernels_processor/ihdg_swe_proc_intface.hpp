#ifndef IHDG_SWE_PROC_INTFACE_HPP
#define IHDG_SWE_PROC_INTFACE_HPP

namespace SWE {
namespace IHDG {
template <typename InterfaceType>
void Problem::init_interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    if (stepper.GetOrder() == 2) {
        const uint stage = stepper.GetStage();

        auto& state_prev_in = intface.data_in.state[stage];
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];

        auto& state_prev_ex = intface.data_ex.state[stage];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            boundary_in.q_at_gp[var] = intface.ComputeUgpIN(state_prev_in.q[var]);
            boundary_ex.q_at_gp[var] = intface.ComputeUgpEX(state_prev_ex.q[var]);
        }
    }
}

template <typename InterfaceType>
void Problem::local_interface_kernel(const ProblemStepperType& stepper, InterfaceType& intface) {
    const uint stage = stepper.GetStage();

    auto& state_in    = intface.data_in.state[stage + 1];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[stage + 1];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        boundary_in.q_at_gp[var] = intface.ComputeUgpIN(state_in.q[var]);
        boundary_ex.q_at_gp[var] = intface.ComputeUgpEX(state_ex.q[var]);
    }
}
}
}

#endif
