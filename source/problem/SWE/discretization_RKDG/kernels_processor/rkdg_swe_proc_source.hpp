#ifndef RKDG_SWE_PROC_SOURCE_HPP
#define RKDG_SWE_PROC_SOURCE_HPP

#include "problem/SWE/problem_function_files/swe_source_functions.hpp"

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::source_kernel(const ProblemStepperType& stepper, ElementType& elt) {
    auto& wd_state = elt.data.wet_dry_state;

    if (wd_state.wet) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;
        // auto& source   = elt.data.source;

        double t = stepper.GetTimeAtCurrentStage();

        // note we assume that the values at gauss points have already been computed
        // compute contribution of hydrostatic pressure
        set_constant(row(internal.source_at_gp, SWE::Variables::ze), 0.0);

        row(internal.source_at_gp, SWE::Variables::qx) =
            Global::g * vec_cw_mult(internal.aux_at_gp[SWE::Auxiliaries::sp],
                                    vec_cw_mult(row(internal.dbath_at_gp, GlobalCoord::x),
                                                internal.q_at_gp[SWE::Variables::ze]));
        row(internal.source_at_gp, SWE::Variables::qy) =
            Global::g *
            vec_cw_mult(row(internal.dbath_at_gp, GlobalCoord::y), internal.q_at_gp[SWE::Variables::ze]);

        if (SWE::SourceTerms::function_source) {
            auto source_q = [t](Point<2>& pt) { return SWE::source_q(t, pt); };

            internal.source_at_gp += elt.ComputeFgp(source_q);
        }

        /*if (SWE::SourceTerms::bottom_friction) {
            double Cf = SWE::SourceTerms::Cf;

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute bottom friction contribution
                double u = internal.q_at_gp(SWE::Variables::qx, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);
                double v = internal.q_at_gp(SWE::Variables::qy, gp) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);

                // compute manning friction factor
                if (source.manning) {
                    Cf = source.g_manning_n_sq / std::pow(internal.aux_at_gp(SWE::Auxiliaries::h, gp), 1.0 / 3.0);
                    if (Cf < SWE::SourceTerms::Cf)
                        Cf = SWE::SourceTerms::Cf;
                }

                double bottom_friction_stress = Cf * std::hypot(u, v) / internal.aux_at_gp(SWE::Auxiliaries::h, gp);

                internal.source_at_gp(SWE::Variables::qx, gp) -=
                    bottom_friction_stress * internal.q_at_gp(SWE::Variables::qx, gp);
                internal.source_at_gp(SWE::Variables::qy, gp) -=
                    bottom_friction_stress * internal.q_at_gp(SWE::Variables::qy, gp);
            }
        }

        if (SWE::SourceTerms::meteo_forcing) {
            // elt.ComputeNodalUgp(source.tau_s, internal.tau_s_at_gp);
            // elt.ComputeNodalDUgp(source.p_atm, internal.dp_atm_at_gp);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute surface friction contribution
                internal.source_at_gp(SWE::Variables::qx, gp) +=
                    internal.tau_s_at_gp(GlobalCoord::x, gp) / Global::rho_water;
                internal.source_at_gp(SWE::Variables::qy, gp) +=
                    internal.tau_s_at_gp(GlobalCoord::y, gp) / Global::rho_water;

                // compute atmospheric pressure contribution
                internal.source_at_gp(SWE::Variables::qx, gp) -=
                    internal.aux_at_gp(SWE::Auxiliaries::sp, gp) * internal.aux_at_gp(SWE::Auxiliaries::h, gp) *
                    internal.dp_atm_at_gp(GlobalCoord::x, gp) / Global::rho_water;
                internal.source_at_gp(SWE::Variables::qy, gp) -= internal.aux_at_gp(SWE::Auxiliaries::h, gp) *
                                                                 internal.dp_atm_at_gp(GlobalCoord::y, gp) /
                                                                 Global::rho_water;
            }
        }

        if (SWE::SourceTerms::tide_potential) {
            // elt.ComputeNodalDUgp(source.tide_pot, internal.dtide_pot_at_gp);

            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute tide potential contribution
                internal.source_at_gp(SWE::Variables::qx, gp) +=
                    Global::g * internal.aux_at_gp(SWE::Auxiliaries::sp, gp) *
                    internal.aux_at_gp(SWE::Auxiliaries::h, gp) * internal.dtide_pot_at_gp(GlobalCoord::x, gp);
                internal.source_at_gp(SWE::Variables::qy, gp) += Global::g *
                                                                 internal.aux_at_gp(SWE::Auxiliaries::h, gp) *
                                                                 internal.dtide_pot_at_gp(GlobalCoord::y, gp);
            }
        }

        if (SWE::SourceTerms::coriolis) {
            for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
                // compute coriolis contribution
                internal.source_at_gp(SWE::Variables::qx, gp) +=
                    source.coriolis_f * internal.q_at_gp(SWE::Variables::qy, gp);
                internal.source_at_gp(SWE::Variables::qy, gp) -=
                    source.coriolis_f * internal.q_at_gp(SWE::Variables::qx, gp);
            }
        }*/
        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            state.rhs[var] += elt.IntegrationPhi(row(internal.source_at_gp,var));
        }
    }
}
}
}

#endif
