#ifndef RKDG_SWE_PROC_VOLUME_HPP
#define RKDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();
    auto& wd_state   = elt.data.wet_dry_state;
    auto& state      = elt.data.state[stage];

    set_constant(state.rhs, 0.0);
    if (wd_state.wet) {
        auto& internal = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);

        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        auto u =
            cwise_division(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
        auto v =
            cwise_division(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

        auto uuh = cwise_multiplication(u, row(internal.q_at_gp, SWE::Variables::qx));
        auto vvh = cwise_multiplication(v, row(internal.q_at_gp, SWE::Variables::qy));
        auto uvh = cwise_multiplication(u, row(internal.q_at_gp, SWE::Variables::qy));
        auto pe  = Global::g * (0.5 * cwise_multiplication(row(internal.q_at_gp, SWE::Variables::ze),
                                                          row(internal.q_at_gp, SWE::Variables::ze)) +
                               cwise_multiplication(row(internal.q_at_gp, SWE::Variables::ze),
                                                    row(internal.aux_at_gp, SWE::Auxiliaries::bath)));

        row(internal.Fx_at_gp, SWE::Variables::ze) = cwise_multiplication(row(internal.aux_at_gp, SWE::Auxiliaries::sp),
                                                                          row(internal.q_at_gp, SWE::Variables::qx));
        row(internal.Fx_at_gp, SWE::Variables::qx) =
            cwise_multiplication(row(internal.aux_at_gp, SWE::Auxiliaries::sp), uuh + pe);
        row(internal.Fx_at_gp, SWE::Variables::qy) =
            cwise_multiplication(row(internal.aux_at_gp, SWE::Auxiliaries::sp), uvh);

        row(internal.Fy_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qy);
        row(internal.Fy_at_gp, SWE::Variables::qx) = uvh;
        row(internal.Fy_at_gp, SWE::Variables::qy) = vvh + pe;

        state.rhs = elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) +
                    elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
    }
}
}
}

#endif
