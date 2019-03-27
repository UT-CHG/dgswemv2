#ifndef RKDG_SWE_PROC_VOLUME_HPP
#define RKDG_SWE_PROC_VOLUME_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace RKDG {
template <typename StepperType, typename ElementType>
void Problem::volume_kernel(const StepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& wd_state = elt.data.wet_dry_state;
    auto& state    = elt.data.state[stage];

    set_constant(state.rhs, 0.0);

    if (wd_state.wet) {
        auto& internal = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);

        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        SWE::get_F(internal.q_at_gp, internal.aux_at_gp, internal.Fx_at_gp, internal.Fy_at_gp);

        // Spherical projection
        row(internal.Fx_at_gp, SWE::Variables::ze) =
            vec_cw_mult(row(internal.aux_at_gp, SWE::Auxiliaries::sp), row(internal.Fx_at_gp, SWE::Variables::ze));
        row(internal.Fx_at_gp, SWE::Variables::qx) =
            vec_cw_mult(row(internal.aux_at_gp, SWE::Auxiliaries::sp), row(internal.Fx_at_gp, SWE::Variables::qx));
        row(internal.Fx_at_gp, SWE::Variables::qy) =
            vec_cw_mult(row(internal.aux_at_gp, SWE::Auxiliaries::sp), row(internal.Fx_at_gp, SWE::Variables::qy));

        state.rhs = elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) +
                    elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
    }
}
}
}

#endif
