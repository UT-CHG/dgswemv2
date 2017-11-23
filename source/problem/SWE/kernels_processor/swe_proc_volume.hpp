#ifndef SWE_PROC_VOLUME_HPP
#define SWE_PROC_VOLUME_HPP

namespace SWE {
template <typename ElementType>
void Problem::volume_kernel(const Stepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.get_stage();

        auto& state = elt.data.state[stage];
        auto& internal = elt.data.internal;

        elt.ComputeUgp(state.ze, internal.ze_at_gp);
        elt.ComputeUgp(state.qx, internal.qx_at_gp);
        elt.ComputeUgp(state.qy, internal.qy_at_gp);

        // assemble flux
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            internal.ze_flux_at_gp[GlobalCoord::x][gp] = internal.qx_at_gp[gp];
            internal.ze_flux_at_gp[GlobalCoord::y][gp] = internal.qy_at_gp[gp];

            internal.qx_flux_at_gp[GlobalCoord::x][gp] = std::pow(internal.qx_at_gp[gp], 2) / internal.h_at_gp[gp] +
                                                         Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) +
                                                                      internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
            internal.qx_flux_at_gp[GlobalCoord::y][gp] =
                internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.h_at_gp[gp];

            internal.qy_flux_at_gp[GlobalCoord::x][gp] =
                internal.qx_at_gp[gp] * internal.qy_at_gp[gp] / internal.h_at_gp[gp];
            internal.qy_flux_at_gp[GlobalCoord::y][gp] = std::pow(internal.qy_at_gp[gp], 2) / internal.h_at_gp[gp] +
                                                         Global::g * (0.5 * std::pow(internal.ze_at_gp[gp], 2) +
                                                                      internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
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
