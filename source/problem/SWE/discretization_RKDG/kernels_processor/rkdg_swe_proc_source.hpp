#ifndef RKDG_SWE_PROC_SOURCE_HPP
#define RKDG_SWE_PROC_SOURCE_HPP

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::source_kernel(const RKStepper& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;
        auto& source   = elt.data.source;
        auto& sp_at_gp = elt.data.spherical_projection.sp_at_gp_internal;

        double t = stepper.GetTimeAtCurrentStage();

        if (SWE::SourceTerms::function_source) {
            auto source_u = [t](Point<2>& pt) { return SWE::source_u(t, pt); };

            elt.ComputeFgp(source_u, internal.source_at_gp);
        } else {
            std::fill(internal.source_at_gp.begin(), internal.source_at_gp.end(), 0.0);
        }

        // note we assume that the values at gauss points have already been computed
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute contribution of hydrostatic pressure
            internal.source_at_gp[gp][SWE::Variables::qx] += Global::g * sp_at_gp[gp] *
                                                             internal.bath_deriv_wrt_x_at_gp[gp] *
                                                             internal.q_at_gp[gp][SWE::Variables::ze];
            internal.source_at_gp[gp][SWE::Variables::qy] +=
                Global::g * internal.bath_deriv_wrt_y_at_gp[gp] * internal.q_at_gp[gp][SWE::Variables::ze];
        }

        if (SWE::SourceTerms::bottom_friction) {
            double Cf = SWE::SourceTerms::Cf;

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute bottom friction contribution
                double u_at_gp = internal.q_at_gp[gp][SWE::Variables::qx] / internal.h_at_gp[gp];
                double v_at_gp = internal.q_at_gp[gp][SWE::Variables::qy] / internal.h_at_gp[gp];

                // compute manning friction factor
                if (source.manning) {
                    Cf = source.g_manning_n_sq / std::pow(internal.h_at_gp[gp], 1.0 / 3.0);
                    if (Cf < SWE::SourceTerms::Cf)
                        Cf = SWE::SourceTerms::Cf;
                }

                double bottom_friction_stress = Cf * std::hypot(u_at_gp, v_at_gp) / internal.h_at_gp[gp];

                internal.source_at_gp[gp][SWE::Variables::qx] -=
                    bottom_friction_stress * internal.q_at_gp[gp][SWE::Variables::qx];
                internal.source_at_gp[gp][SWE::Variables::qy] -=
                    bottom_friction_stress * internal.q_at_gp[gp][SWE::Variables::qy];
            }
        }

        if (SWE::SourceTerms::meteo_forcing) {
            elt.ComputeNodalUgp(source.tau_s[GlobalCoord::x], internal.tau_s_at_gp[GlobalCoord::x]);
            elt.ComputeNodalUgp(source.tau_s[GlobalCoord::y], internal.tau_s_at_gp[GlobalCoord::y]);

            elt.ComputeNodalDUgp(GlobalCoord::x, source.p_atm, internal.dp_atm_at_gp[GlobalCoord::x]);
            elt.ComputeNodalDUgp(GlobalCoord::y, source.p_atm, internal.dp_atm_at_gp[GlobalCoord::y]);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute surface friction contribution
                internal.source_at_gp[gp][SWE::Variables::qx] +=
                    internal.tau_s_at_gp[GlobalCoord::x][gp] / Global::rho_water;
                internal.source_at_gp[gp][SWE::Variables::qy] +=
                    internal.tau_s_at_gp[GlobalCoord::y][gp] / Global::rho_water;

                // compute atmospheric pressure contribution
                internal.source_at_gp[gp][SWE::Variables::qx] -=
                    sp_at_gp[gp] * internal.h_at_gp[gp] * internal.dp_atm_at_gp[GlobalCoord::x][gp] / Global::rho_water;
                internal.source_at_gp[gp][SWE::Variables::qy] -=
                    internal.h_at_gp[gp] * internal.dp_atm_at_gp[GlobalCoord::y][gp] / Global::rho_water;
            }
        }

        if (SWE::SourceTerms::tide_potential) {
            elt.ComputeNodalDUgp(GlobalCoord::x, source.tide_pot, internal.dtide_pot_at_gp[GlobalCoord::x]);
            elt.ComputeNodalDUgp(GlobalCoord::y, source.tide_pot, internal.dtide_pot_at_gp[GlobalCoord::y]);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute tide potential contribution
                internal.source_at_gp[gp][SWE::Variables::qx] +=
                    Global::g * sp_at_gp[gp] * internal.h_at_gp[gp] * internal.dtide_pot_at_gp[GlobalCoord::x][gp];
                internal.source_at_gp[gp][SWE::Variables::qy] +=
                    Global::g * internal.h_at_gp[gp] * internal.dtide_pot_at_gp[GlobalCoord::y][gp];
            }
        }

        if (SWE::SourceTerms::coriolis) {
            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute coriolis contribution
                internal.source_at_gp[gp][SWE::Variables::qx] +=
                    source.coriolis_f * internal.q_at_gp[gp][SWE::Variables::qy];
                internal.source_at_gp[gp][SWE::Variables::qy] -=
                    source.coriolis_f * internal.q_at_gp[gp][SWE::Variables::qx];
            }
        }

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.rhs[dof] += elt.IntegrationPhi(dof, internal.source_at_gp);
        }
    }
}
}
}

#endif
