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
        auto& sp_at_gp = elt.data.spherical_projection.sp_at_gp_internal;

        elt.ComputeUgp(state.q, internal.q_at_gp);

        double u_at_gp = 0.0;
        double v_at_gp = 0.0;

        double uuh_at_gp = 0.0;
        double vvh_at_gp = 0.0;
        double uvh_at_gp = 0.0;
        double pe_at_gp  = 0.0;

        // assemble flux
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            internal.h_at_gp[gp] = internal.q_at_gp[gp][SWE::Variables::ze] + internal.bath_at_gp[gp];

            u_at_gp = internal.q_at_gp[gp][SWE::Variables::qx] / internal.h_at_gp[gp];
            v_at_gp = internal.q_at_gp[gp][SWE::Variables::qy] / internal.h_at_gp[gp];

            uuh_at_gp = u_at_gp * internal.q_at_gp[gp][SWE::Variables::qx];
            vvh_at_gp = v_at_gp * internal.q_at_gp[gp][SWE::Variables::qy];
            uvh_at_gp = u_at_gp * internal.q_at_gp[gp][SWE::Variables::qy];
            pe_at_gp  = Global::g * (0.5 * std::pow(internal.q_at_gp[gp][SWE::Variables::ze], 2) +
                                    internal.q_at_gp[gp][SWE::Variables::ze] * internal.bath_at_gp[gp]);

            internal.Fx_at_gp[gp][SWE::Variables::ze] = sp_at_gp[gp] * internal.q_at_gp[gp][SWE::Variables::qx];
            internal.Fx_at_gp[gp][SWE::Variables::qx] = sp_at_gp[gp] * (uuh_at_gp + pe_at_gp);
            internal.Fx_at_gp[gp][SWE::Variables::qy] = sp_at_gp[gp] * uvh_at_gp;

            internal.Fy_at_gp[gp][SWE::Variables::ze] = internal.q_at_gp[gp][SWE::Variables::qy];
            internal.Fy_at_gp[gp][SWE::Variables::qx] = uvh_at_gp;
            internal.Fy_at_gp[gp][SWE::Variables::qy] = vvh_at_gp + pe_at_gp;
        }

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.rhs[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.Fx_at_gp) +
                             elt.IntegrationDPhi(GlobalCoord::y, dof, internal.Fy_at_gp);
        }
    }
}
}
}

#endif
