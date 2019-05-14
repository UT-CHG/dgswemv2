#ifndef EHDG_SWE_PROC_VOLUME_HPP
#define EHDG_SWE_PROC_VOLUME_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace EHDG {
template <typename ElementType>
void Problem::local_volume_kernel(const ProblemStepperType& stepper, ElementType& elt) {
    auto& state = elt.data.state[stepper.GetStage()];

    for ( uint var = 0; var < SWE::n_variables; ++ var ) {
        set_constant(state.rhs[var], 0.0);
    }

    if (elt.data.wet_dry_state.wet) {
        auto& internal = elt.data.internal;

        internal.q_at_gp[SWE::Variables::ze] = elt.ComputeUgp(state.q[SWE::Variables::ze]);
        internal.q_at_gp[SWE::Variables::qx] = elt.ComputeUgp(state.q[SWE::Variables::qx]);
        internal.q_at_gp[SWE::Variables::qy] = elt.ComputeUgp(state.q[SWE::Variables::qy]);

        internal.aux_at_gp[SWE::Auxiliaries::h] =
            internal.q_at_gp[SWE::Variables::ze] + internal.aux_at_gp[SWE::Auxiliaries::bath];

        SWE::get_F(internal.q_at_gp, internal.aux_at_gp, internal.Fx_at_gp, internal.Fy_at_gp);

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            state.rhs[var] = elt.IntegrationDPhi(GlobalCoord::x, row(internal.Fx_at_gp, var)) +
                elt.IntegrationDPhi(GlobalCoord::y, row(internal.Fy_at_gp, var));
        }
    }
}
}
}

#endif
