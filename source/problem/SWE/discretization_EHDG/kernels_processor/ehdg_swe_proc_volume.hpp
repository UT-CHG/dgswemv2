#ifndef EHDG_SWE_PROC_VOLUME_HPP
#define EHDG_SWE_PROC_VOLUME_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace EHDG {
template <typename StepperType, typename ElementType>
void Problem::local_volume_kernel(const StepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    row(internal.aux_at_gp, SWE::Auxiliaries::h) =
        row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

    SWE::get_F(internal.q_at_gp, internal.aux_at_gp, internal.Fx_at_gp, internal.Fy_at_gp);

    state.rhs =
        elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) + elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
}
}
}

#endif
