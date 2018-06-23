#ifndef EHDG_SWE_PROC_VOLUME_HPP
#define EHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace EHDG {
template <typename ElementType>
void Problem::local_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    elt.ComputeUgp(state.ze, internal.ze_at_gp);
    elt.ComputeUgp(state.qx, internal.qx_at_gp);
    elt.ComputeUgp(state.qy, internal.qy_at_gp);

    double u_at_gp = 0.0;
    double v_at_gp = 0.0;

    double uuh_at_gp = 0.0;
    double vvh_at_gp = 0.0;
    double uvh_at_gp = 0.0;
    double pe_at_gp  = 0.0;

    // assemble flux
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.h_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];

        u_at_gp = internal.qx_at_gp[gp] / internal.h_at_gp[gp];
        v_at_gp = internal.qy_at_gp[gp] / internal.h_at_gp[gp];

        uuh_at_gp = u_at_gp * internal.qx_at_gp[gp];
        vvh_at_gp = v_at_gp * internal.qy_at_gp[gp];
        uvh_at_gp = u_at_gp * internal.qy_at_gp[gp];
        pe_at_gp =
            Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);

        internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
        internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];

        internal.qx_flux_at_gp[GlobalCoord::x][gp] = (uuh_at_gp + pe_at_gp);
        internal.qx_flux_at_gp[GlobalCoord::y][gp] = uvh_at_gp;

        internal.qy_flux_at_gp[GlobalCoord::x][gp] = uvh_at_gp;
        internal.qy_flux_at_gp[GlobalCoord::y][gp] = vvh_at_gp + pe_at_gp;
    }

    // skip dof = 0, which is a constant and thus trivially 0 NOT ALWAYS!
    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        state.rhs_ze[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.ze_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.ze_flux_at_gp[GlobalCoord::y]);

        state.rhs_qx[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qx_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qx_flux_at_gp[GlobalCoord::y]);

        state.rhs_qy[dof] = elt.IntegrationDPhi(GlobalCoord::x, dof, internal.qy_flux_at_gp[GlobalCoord::x]) +
                            elt.IntegrationDPhi(GlobalCoord::y, dof, internal.qy_flux_at_gp[GlobalCoord::y]);
    }
}
}
}

#endif
