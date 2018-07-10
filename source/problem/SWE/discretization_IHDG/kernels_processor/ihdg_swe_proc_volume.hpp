#ifndef IHDG_SWE_PROC_VOLUME_HPP
#define IHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace IHDG {
template <typename ElementType>
void Problem::prepare_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state_prev = elt.data.state[stage];
    auto& state      = elt.data.state[stage + 1];
    auto& internal   = elt.data.internal;

    elt.ComputeUgp(state.q, internal.q_at_gp);
    elt.ComputeUgp(state_prev.q, internal.q_prev_at_gp);

    double u_at_gp = 0.0;
    double v_at_gp = 0.0;

    double uuh_at_gp = 0.0;
    double vvh_at_gp = 0.0;
    double uvh_at_gp = 0.0;
    double pe_at_gp  = 0.0;

    // assemble flux
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.aux_at_gp[gp][SWE::Auxiliaries::h] =
            internal.q_at_gp[gp][SWE::Variables::ze] + internal.aux_at_gp[gp][SWE::Auxiliaries::bath];

        u_at_gp = internal.q_at_gp[gp][SWE::Variables::qx] / internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        v_at_gp = internal.q_at_gp[gp][SWE::Variables::qy] / internal.aux_at_gp[gp][SWE::Auxiliaries::h];

        uuh_at_gp = u_at_gp * internal.q_at_gp[gp][SWE::Variables::qx];
        vvh_at_gp = v_at_gp * internal.q_at_gp[gp][SWE::Variables::qy];
        uvh_at_gp = u_at_gp * internal.q_at_gp[gp][SWE::Variables::qy];
        pe_at_gp =
            Global::g * (0.5 * std::pow(internal.q_at_gp[gp][SWE::Variables::ze], 2) +
                         internal.q_at_gp[gp][SWE::Variables::ze] * internal.aux_at_gp[gp][SWE::Auxiliaries::bath]);

        // Flux terms
        internal.Fx_at_gp[gp][SWE::Variables::ze] = internal.q_at_gp[gp][SWE::Variables::qx];
        internal.Fx_at_gp[gp][SWE::Variables::qx] = (uuh_at_gp + pe_at_gp);
        internal.Fx_at_gp[gp][SWE::Variables::qy] = uvh_at_gp;

        internal.Fy_at_gp[gp][SWE::Variables::ze] = internal.q_at_gp[gp][SWE::Variables::qy];
        internal.Fy_at_gp[gp][SWE::Variables::qx] = uvh_at_gp;
        internal.Fy_at_gp[gp][SWE::Variables::qy] = vvh_at_gp + pe_at_gp;

        // dFx/dq terms
        internal.dFx_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        internal.dFx_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qx) = 1.0;
        internal.dFx_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qy) = 0.0;

        internal.dFx_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::ze) =
            -u_at_gp * u_at_gp + Global::g * internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        internal.dFx_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qx) = 2 * u_at_gp;
        internal.dFx_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qy) = 0.0;

        internal.dFx_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::ze) = -u_at_gp * v_at_gp;
        internal.dFx_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qx) = v_at_gp;
        internal.dFx_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qy) = u_at_gp;

        // dFy/dq terms
        internal.dFy_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::ze) = 0.0;
        internal.dFy_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qx) = 0.0;
        internal.dFy_dq_at_gp[gp](SWE::Variables::ze, SWE::Variables::qy) = 1.0;

        internal.dFy_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::ze) = -u_at_gp * v_at_gp;
        internal.dFy_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qx) = v_at_gp;
        internal.dFy_dq_at_gp[gp](SWE::Variables::qx, SWE::Variables::qy) = u_at_gp;

        internal.dFy_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::ze) =
            -v_at_gp * v_at_gp + Global::g * internal.aux_at_gp[gp][SWE::Auxiliaries::h];
        internal.dFy_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qx) = 0.0;
        internal.dFy_dq_at_gp[gp](SWE::Variables::qy, SWE::Variables::qy) = 2 * v_at_gp;
    }
}
}
}

#endif
