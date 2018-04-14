#ifndef SWE_PROC_SOURCE_HPP
#define SWE_PROC_SOURCE_HPP

namespace SWE {
template <typename ElementType>
void Problem::source_kernel(const Stepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;
        auto& source   = elt.data.source;
        auto& sp_at_gp = elt.data.spherical_projection.sp_at_gp_internal;

        double t = stepper.GetTimeAtCurrentStage();

        if (SWE::SourceTerms::function_source) {
            auto source_ze = [t](Point<2>& pt) { return SWE::source_ze(t, pt); };

            auto source_qx = [t](Point<2>& pt) { return SWE::source_qx(t, pt); };

            auto source_qy = [t](Point<2>& pt) { return SWE::source_qy(t, pt); };

            elt.ComputeFgp(source_ze, internal.ze_source_term_at_gp);
            elt.ComputeFgp(source_qx, internal.qx_source_term_at_gp);
            elt.ComputeFgp(source_qy, internal.qy_source_term_at_gp);
        } else {
            std::fill(internal.ze_source_term_at_gp.begin(), internal.ze_source_term_at_gp.end(), 0.0);
            std::fill(internal.qx_source_term_at_gp.begin(), internal.qx_source_term_at_gp.end(), 0.0);
            std::fill(internal.qy_source_term_at_gp.begin(), internal.qy_source_term_at_gp.end(), 0.0);
        }

        // note we assume that the values at gauss points have already been computed
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute contribution of hydrostatic pressure
            internal.qx_source_term_at_gp[gp] +=
                Global::g * sp_at_gp[gp] * internal.bath_deriv_wrt_x_at_gp[gp] * internal.ze_at_gp[gp];
            internal.qy_source_term_at_gp[gp] +=
                Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.ze_at_gp[gp];
        }

        if (SWE::SourceTerms::bottom_friction) {
            double Cf = Global::Cf;

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute bottom friction contribution
                double u_at_gp = internal.qx_at_gp[gp] / internal.h_at_gp[gp];
                double v_at_gp = internal.qy_at_gp[gp] / internal.h_at_gp[gp];

                // compute manning friction factor
                if (source.manning) {
                    Cf = source.g_manning_n_sq / std::pow(internal.h_at_gp[gp], 1.0 / 3.0);
                    if (Cf < Global::Cf)
                        Cf = Global::Cf;
                }

                double bottom_friction_stress = Cf * std::hypot(u_at_gp, v_at_gp) / internal.h_at_gp[gp];

                internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
                internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
            }
        }

        if (SWE::SourceTerms::meteo_forcing) {
            elt.ComputeNodalUgp(source.tau_s[GlobalCoord::x], internal.tau_s_at_gp[GlobalCoord::x]);
            elt.ComputeNodalUgp(source.tau_s[GlobalCoord::y], internal.tau_s_at_gp[GlobalCoord::y]);

            elt.ComputeNodalDUgp(GlobalCoord::x, source.p_atm, internal.dp_atm_at_gp[GlobalCoord::x]);
            elt.ComputeNodalDUgp(GlobalCoord::y, source.p_atm, internal.dp_atm_at_gp[GlobalCoord::y]);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute surface friction contribution
                internal.qx_source_term_at_gp[gp] += internal.tau_s_at_gp[GlobalCoord::x][gp] / Global::rho_water;
                internal.qy_source_term_at_gp[gp] += internal.tau_s_at_gp[GlobalCoord::y][gp] / Global::rho_water;

                // compute atmospheric pressure contribution
                internal.qx_source_term_at_gp[gp] -=
                    sp_at_gp[gp] * internal.h_at_gp[gp] * internal.dp_atm_at_gp[GlobalCoord::x][gp] / Global::rho_water;
                internal.qy_source_term_at_gp[gp] -=
                    internal.h_at_gp[gp] * internal.dp_atm_at_gp[GlobalCoord::y][gp] / Global::rho_water;
            }
        }

        if (SWE::SourceTerms::tidal_potential) {
            elt.ComputeNodalDUgp(GlobalCoord::x, source.tidal_pot, internal.dtidal_pot_at_gp[GlobalCoord::x]);
            elt.ComputeNodalDUgp(GlobalCoord::y, source.tidal_pot, internal.dtidal_pot_at_gp[GlobalCoord::y]);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute tidal potential contribution
                internal.qx_source_term_at_gp[gp] +=
                    Global::g * sp_at_gp[gp] * internal.h_at_gp[gp] * internal.dtidal_pot_at_gp[GlobalCoord::x][gp];
                internal.qy_source_term_at_gp[gp] +=
                    Global::g * internal.h_at_gp[gp] * internal.dtidal_pot_at_gp[GlobalCoord::y][gp];
            }
        }

        if (SWE::SourceTerms::coriolis) {
            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute coriolis contribution
                internal.qx_source_term_at_gp[gp] += source.coriolis_f * internal.qy_at_gp[gp];
                internal.qy_source_term_at_gp[gp] -= source.coriolis_f * internal.qx_at_gp[gp];
            }
        }

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.rhs_ze[dof] += elt.IntegrationPhi(dof, internal.ze_source_term_at_gp);
            state.rhs_qx[dof] += elt.IntegrationPhi(dof, internal.qx_source_term_at_gp);
            state.rhs_qy[dof] += elt.IntegrationPhi(dof, internal.qy_source_term_at_gp);
        }
    }
}
}

#endif
