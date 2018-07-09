#ifndef EHDG_SWE_PROC_VOLUME_HPP
#define EHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace EHDG {
template <typename ElementType>
void Problem::local_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    elt.ComputeUgp(state.q, internal.q_at_gp);

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

        internal.Fx_at_gp[gp][SWE::Variables::ze] = internal.q_at_gp[gp][SWE::Variables::qx];
        internal.Fx_at_gp[gp][SWE::Variables::qx] = (uuh_at_gp + pe_at_gp);
        internal.Fx_at_gp[gp][SWE::Variables::qy] = uvh_at_gp;

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

#endif
