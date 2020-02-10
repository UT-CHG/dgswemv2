#ifndef SWE_SOURCE_HPP
#define SWE_SOURCE_HPP

#include "problem/SWE/problem_data_structure/swe_data_source.hpp"
#include "problem/SWE/problem_function_files/swe_source_functions.hpp"
#include "problem/SWE/seabed_update/swe_seabed_update.hpp"

namespace SWE {
template <typename StepperType, typename ElementType>
void get_source(const StepperType& stepper, ElementType& elt) {
    auto& state    = elt.data.state[stepper.GetStage()];
    auto& internal = elt.data.internal;
    auto& source   = elt.data.source;

    // note we assume that the values at gauss points have already been computed
    // compute contribution of hydrostatic pressure
    row(internal.dh_at_gp, GlobalCoord::x) =
        elt.ComputeDUgp(GlobalCoord::x, row(state.q, SWE::Variables::ze)) + row(internal.db_at_gp, GlobalCoord::x);
    row(internal.dh_at_gp, GlobalCoord::y) =
        elt.ComputeDUgp(GlobalCoord::y, row(state.q, SWE::Variables::ze)) + row(internal.db_at_gp, GlobalCoord::y);

    row(internal.dhc_at_gp, GlobalCoord::x) = elt.ComputeDUgp(GlobalCoord::x, row(state.q, SWE::Variables::hc));
    row(internal.dhc_at_gp, GlobalCoord::y) = elt.ComputeDUgp(GlobalCoord::y, row(state.q, SWE::Variables::hc));

    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        internal.rho_mix_at_gp[gp] = rho_mixture(column(internal.q_at_gp, gp), column(internal.aux_at_gp, gp));
        internal.E_at_gp[gp]       = entrainment_rate(column(internal.q_at_gp, gp), column(internal.aux_at_gp, gp));
        internal.D_at_gp[gp]       = deposition_rate(column(internal.q_at_gp, gp), column(internal.aux_at_gp, gp));
        internal.sed_mom_frac_at_gp[gp] = (Global::rho_bed - internal.rho_mix_at_gp[gp]) /
                                          (internal.rho_mix_at_gp[gp] * (1.0 - Global::sat_sediment));
    }

    const auto b_grad_x =
        Global::g * vec_cw_mult(row(internal.db_at_gp, GlobalCoord::x), row(internal.q_at_gp, SWE::Variables::ze));
    const auto b_grad_y =
        Global::g * vec_cw_mult(row(internal.db_at_gp, GlobalCoord::y), row(internal.q_at_gp, SWE::Variables::ze));

    const auto c_grad_x =
        -Global::g * (Global::rho_sediment - Global::rho_water) / 2.0 *
        vec_cw_div(vec_cw_mult(row(internal.dhc_at_gp, GlobalCoord::x), row(internal.aux_at_gp, SWE::Auxiliaries::h)) -
                       vec_cw_mult(row(internal.dh_at_gp, GlobalCoord::x), row(internal.q_at_gp, SWE::Variables::hc)),
                   internal.rho_mix_at_gp);
    const auto c_grad_y =
        -Global::g * (Global::rho_sediment - Global::rho_water) / 2.0 *
        vec_cw_div(vec_cw_mult(row(internal.dhc_at_gp, GlobalCoord::y), row(internal.aux_at_gp, SWE::Auxiliaries::h)) -
                       vec_cw_mult(row(internal.dh_at_gp, GlobalCoord::y), row(internal.q_at_gp, SWE::Variables::hc)),
                   internal.rho_mix_at_gp);

    const auto u = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
    const auto v = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));
    const auto sed_mom_x =
        -vec_cw_mult(vec_cw_mult(internal.sed_mom_frac_at_gp, internal.E_at_gp - internal.D_at_gp), u);
    const auto sed_mom_y =
        -vec_cw_mult(vec_cw_mult(internal.sed_mom_frac_at_gp, internal.E_at_gp - internal.D_at_gp), v);

    row(internal.source_at_gp, SWE::Variables::ze) =
        1.0 / (1.0 - Global::sat_sediment) * (internal.E_at_gp - internal.D_at_gp);
    row(internal.source_at_gp, SWE::Variables::qx) = b_grad_x + c_grad_x + sed_mom_x;
    row(internal.source_at_gp, SWE::Variables::qy) = b_grad_y + c_grad_y + sed_mom_y;
    row(internal.source_at_gp, SWE::Variables::hc) = internal.E_at_gp - internal.D_at_gp;

    if (SWE::SourceTerms::function_source) {
        const double t = stepper.GetTimeAtCurrentStage();
        internal.source_at_gp += elt.ComputeFgp([t](Point<2>& pt) { return SWE::source_q(t, pt); });
    }

    if (SWE::SourceTerms::bottom_friction) {
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

    /*if (SWE::SourceTerms::meteo_forcing) {
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
}
}

#endif