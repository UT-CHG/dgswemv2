#ifndef RKDG_SWE_PROC_VOLUME_HPP
#define RKDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::volume_kernel(const RKStepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);

        double u = 0.0;
        double v = 0.0;

        double uuh = 0.0;
        double vvh = 0.0;
        double uvh = 0.0;
        double pe  = 0.0;

        // assemble flux
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            internal.aux_at_gp(SWE::Auxiliaries::h, gp) =
                internal.q_at_gp(SWE::Variables::ze, gp) + internal.aux_at_gp(SWE::Auxiliaries::bath, gp);

            u = internal.q_at_gp(SWE::Variables::qx, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);
            v = internal.q_at_gp(SWE::Variables::qy, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);

            uuh = u * internal.q_at_gp(SWE::Variables::qx, gp);
            vvh = v * internal.q_at_gp(SWE::Variables::qy, gp);
            uvh = u * internal.q_at_gp(SWE::Variables::qy, gp);
            pe =
                Global::g * (0.5 * std::pow(internal.q_at_gp(SWE::Variables::ze, gp), 2) +
                             internal.q_at_gp(SWE::Variables::ze, gp) * internal.aux_at_gp(SWE::Auxiliaries::bath, gp));

            internal.Fx_at_gp(SWE::Variables::ze, gp) =
                internal.aux_at_gp(SWE::Auxiliaries::sp, gp) * internal.q_at_gp(SWE::Variables::qx, gp);
            internal.Fx_at_gp(SWE::Variables::qx, gp) = internal.aux_at_gp(SWE::Auxiliaries::sp, gp) * (uuh + pe);
            internal.Fx_at_gp(SWE::Variables::qy, gp) = internal.aux_at_gp(SWE::Auxiliaries::sp, gp) * uvh;

            internal.Fy_at_gp(SWE::Variables::ze, gp) = internal.q_at_gp(SWE::Variables::qy, gp);
            internal.Fy_at_gp(SWE::Variables::qx, gp) = uvh;
            internal.Fy_at_gp(SWE::Variables::qy, gp) = vvh + pe;
        }

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            column(state.rhs, dof) = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.Fx_at_gp) +
                                     elt.IntegrationDPhi(GlobalCoord::y, dof, internal.Fy_at_gp);
        }
    }
}
}
}

#endif
