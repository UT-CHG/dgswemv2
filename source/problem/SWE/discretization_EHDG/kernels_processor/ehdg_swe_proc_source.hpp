#ifndef EHDG_SWE_PROC_SOURCE_HPP
#define EHDG_SWE_PROC_SOURCE_HPP

namespace SWE {
namespace EHDG {
template <typename ElementType>
void Problem::local_source_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;
    auto& source   = elt.data.source;

    double t = stepper.GetTimeAtCurrentStage();

    if (SWE::SourceTerms::function_source) {
        auto source_u = [t](Point<2>& pt) { return SWE::source_u(t, pt); };

        internal.source_at_gp = elt.ComputeFgp(source_u);
    } else {
        internal.source_at_gp = 0.0;
    }

    // note we assume that the values at gauss points have already been computed
    for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
        // compute contribution of hydrostatic pressure
        internal.source_at_gp(SWE::Variables::qx, gp) +=
            Global::g * internal.dbath_at_gp(GlobalCoord::x, gp) * internal.q_at_gp(SWE::Variables::ze, gp);
        internal.source_at_gp(SWE::Variables::qy, gp) +=
            Global::g * internal.dbath_at_gp(GlobalCoord::y, gp) * internal.q_at_gp(SWE::Variables::ze, gp);
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

    if (SWE::SourceTerms::meteo_forcing) {
        // internal.tau_s_at_gp = elt.ComputeNodalUgp(source.tau_s);
        // row(internal.dp_atm_at_gp, GlobalCoord::x) = elt.ComputeNodalDUgp(GlobalCoord::x, source.p_atm);
        // row(internal.dp_atm_at_gp, GlobalCoord::y) = elt.ComputeNodalDUgp(GlobalCoord::y, source.p_atm);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute surface friction contribution
            internal.source_at_gp(SWE::Variables::qx, gp) +=
                internal.tau_s_at_gp(GlobalCoord::x, gp) / Global::rho_water;
            internal.source_at_gp(SWE::Variables::qy, gp) +=
                internal.tau_s_at_gp(GlobalCoord::y, gp) / Global::rho_water;

            // compute atmospheric pressure contribution
            internal.source_at_gp(SWE::Variables::qx, gp) -= internal.aux_at_gp(SWE::Auxiliaries::h, gp) *
                                                             internal.dp_atm_at_gp(GlobalCoord::x, gp) /
                                                             Global::rho_water;
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
                Global::g * internal.aux_at_gp(SWE::Auxiliaries::h, gp) * internal.dtide_pot_at_gp(GlobalCoord::x, gp);
            internal.source_at_gp(SWE::Variables::qy, gp) +=
                Global::g * internal.aux_at_gp(SWE::Auxiliaries::h, gp) * internal.dtide_pot_at_gp(GlobalCoord::y, gp);
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
    }

    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        column(state.rhs, dof) += elt.IntegrationPhi(dof, internal.source_at_gp);
    }
}
}
}

#endif
