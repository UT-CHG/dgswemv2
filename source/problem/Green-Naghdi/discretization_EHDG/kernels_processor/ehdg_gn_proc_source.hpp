#ifndef EHDG_GN_PROC_SOURCE_HPP
#define EHDG_GN_PROC_SOURCE_HPP

namespace GN {
namespace EHDG {
template <typename ElementType>
void Problem::local_swe_source_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;
    // auto& source   = elt.data.source;

    double t = stepper.GetTimeAtCurrentStage();

    // note we assume that the values at gauss points have already been computed
    // compute contribution of hydrostatic pressure
    set_constant(row(internal.source_at_gp, GN::Variables::ze), 0.0);

    row(internal.source_at_gp, GN::Variables::qx) =
        Global::g *
        cwise_multiplication(row(internal.dbath_at_gp, GlobalCoord::x), row(internal.q_at_gp, GN::Variables::ze));

    row(internal.source_at_gp, GN::Variables::qy) =
        Global::g *
        cwise_multiplication(row(internal.dbath_at_gp, GlobalCoord::y), row(internal.q_at_gp, GN::Variables::ze));

    if (GN::SourceTerms::function_source) {
        auto source_u = [t](Point<2>& pt) { return GN::source_u(t, pt); };

        internal.source_at_gp += elt.ComputeFgp(source_u);
    }

    /*if (GN::SourceTerms::bottom_friction) {
        double Cf = GN::SourceTerms::Cf;

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute bottom friction contribution
            double u = internal.q_at_gp(GN::Variables::qx, gp) / internal.aux_at_gp(GN::Auxiliaries::h, gp);
            double v = internal.q_at_gp(GN::Variables::qy, gp) / internal.aux_at_gp(GN::Auxiliaries::h, gp);

            // compute manning friction factor
            if (source.manning) {
                Cf = source.g_manning_n_sq / std::pow(internal.aux_at_gp(GN::Auxiliaries::h, gp), 1.0 / 3.0);
                if (Cf < GN::SourceTerms::Cf)
                    Cf = GN::SourceTerms::Cf;
            }

            double bottom_friction_stress = Cf * std::hypot(u, v) / internal.aux_at_gp(GN::Auxiliaries::h, gp);

            internal.source_at_gp(GN::Variables::qx, gp) -=
                bottom_friction_stress * internal.q_at_gp(GN::Variables::qx, gp);
            internal.source_at_gp(GN::Variables::qy, gp) -=
                bottom_friction_stress * internal.q_at_gp(GN::Variables::qy, gp);
        }
    }

    if (GN::SourceTerms::meteo_forcing) {
        // elt.ComputeNodalUgp(source.tau_s, internal.tau_s_at_gp);
        // elt.ComputeNodalDUgp(source.p_atm, internal.dp_atm_at_gp);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute surface friction contribution
            internal.source_at_gp(GN::Variables::qx, gp) +=
                internal.tau_s_at_gp(GlobalCoord::x, gp) / Global::rho_water;
            internal.source_at_gp(GN::Variables::qy, gp) +=
                internal.tau_s_at_gp(GlobalCoord::y, gp) / Global::rho_water;

            // compute atmospheric pressure contribution
            internal.source_at_gp(GN::Variables::qx, gp) -=
                internal.aux_at_gp(GN::Auxiliaries::sp, gp) * internal.aux_at_gp(GN::Auxiliaries::h, gp) *
                internal.dp_atm_at_gp(GlobalCoord::x, gp) / Global::rho_water;
            internal.source_at_gp(GN::Variables::qy, gp) -= internal.aux_at_gp(GN::Auxiliaries::h, gp) *
                                                             internal.dp_atm_at_gp(GlobalCoord::y, gp) /
                                                             Global::rho_water;
        }
    }

    if (GN::SourceTerms::tide_potential) {
        // elt.ComputeNodalDUgp(source.tide_pot, internal.dtide_pot_at_gp);

        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute tide potential contribution
            internal.source_at_gp(GN::Variables::qx, gp) +=
                Global::g * internal.aux_at_gp(GN::Auxiliaries::sp, gp) *
                internal.aux_at_gp(GN::Auxiliaries::h, gp) * internal.dtide_pot_at_gp(GlobalCoord::x, gp);
            internal.source_at_gp(GN::Variables::qy, gp) += Global::g *
                                                             internal.aux_at_gp(GN::Auxiliaries::h, gp) *
                                                             internal.dtide_pot_at_gp(GlobalCoord::y, gp);
        }
    }

    if (GN::SourceTerms::coriolis) {
        for (uint gp = 0; gp < elt.data.get_ngp_internal(); ++gp) {
            // compute coriolis contribution
            internal.source_at_gp(GN::Variables::qx, gp) +=
                source.coriolis_f * internal.q_at_gp(GN::Variables::qy, gp);
            internal.source_at_gp(GN::Variables::qy, gp) -=
                source.coriolis_f * internal.q_at_gp(GN::Variables::qx, gp);
        }
    }*/

    state.rhs += elt.IntegrationPhi(internal.source_at_gp);
}

template <typename ElementType>
void Problem::dc_source_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    // at this point h_at_gp, u_at_gp and du_at_gp
    // have been calculated in derivatives kernel
    internal.ddu_at_gp = elt.ComputeUgp(state.ddu);

    auto dze_dx = elt.ComputeDUgp(GlobalCoord::x, row(state.q, GN::Variables::ze));
    auto dze_dy = elt.ComputeDUgp(GlobalCoord::y, row(state.q, GN::Variables::ze));

    auto dh_dx = dze_dx + row(internal.dbath_at_gp, GlobalCoord::x);
    auto dh_dy = dze_dy + row(internal.dbath_at_gp, GlobalCoord::y);

    auto h  = row(internal.aux_at_gp, GN::Auxiliaries::h);
    auto h2 = cwise_multiplication(h, h);
    auto h3 = cwise_multiplication(h2, h);

    auto ux = row(internal.du_at_gp, GN::FirstDerivatives::ux);
    auto uy = row(internal.du_at_gp, GN::FirstDerivatives::uy);
    auto vx = row(internal.du_at_gp, GN::FirstDerivatives::vx);
    auto vy = row(internal.du_at_gp, GN::FirstDerivatives::vy);

    auto uxx = row(internal.ddu_at_gp, GN::SecondDerivatives::uxx);
    auto uxy = row(internal.ddu_at_gp, GN::SecondDerivatives::uxy);
    auto uyx = row(internal.ddu_at_gp, GN::SecondDerivatives::uyx);
    auto uyy = row(internal.ddu_at_gp, GN::SecondDerivatives::uyy);
    auto vxx = row(internal.ddu_at_gp, GN::SecondDerivatives::vxx);
    auto vxy = row(internal.ddu_at_gp, GN::SecondDerivatives::vxy);
    auto vyx = row(internal.ddu_at_gp, GN::SecondDerivatives::vyx);
    auto vyy = row(internal.ddu_at_gp, GN::SecondDerivatives::vyy);

    auto c1 = cwise_multiplication(vx, uy) + cwise_multiplication(ux, ux) + cwise_multiplication(ux, vy) +
              cwise_multiplication(vy, vy);

    auto c2 = -cwise_multiplication(vy, uxx) - cwise_multiplication(ux, vyx) + cwise_multiplication(uy, vxx) +
              cwise_multiplication(vx, uyx) + 2.0 * cwise_multiplication(ux, uxx) +
              2.0 * cwise_multiplication(vy, uxx) + 2.0 * cwise_multiplication(ux, vyx) +
              2.0 * cwise_multiplication(vy, vyx);

    auto c3 = -cwise_multiplication(vy, uxy) - cwise_multiplication(ux, vyy) + cwise_multiplication(uy, vxy) +
              cwise_multiplication(vx, uyy) + 2.0 * cwise_multiplication(ux, uxy) +
              2.0 * cwise_multiplication(vy, uxy) + 2.0 * cwise_multiplication(ux, vyy) +
              2.0 * cwise_multiplication(vy, vyy);

    row(internal.w1_rhs_kernel_at_gp, GlobalCoord::x) =
        2.0 * cwise_multiplication(cwise_multiplication(h2, c1), dh_dx) + 2.0 / 3.0 * cwise_multiplication(h3, c2) +
        Global::g / NDParameters::alpha * cwise_multiplication(dze_dx, h);

    row(internal.w1_rhs_kernel_at_gp, GlobalCoord::y) =
        2.0 * cwise_multiplication(cwise_multiplication(h2, c1), dh_dy) + 2.0 / 3.0 * cwise_multiplication(h3, c3) +
        Global::g / NDParameters::alpha * cwise_multiplication(dze_dy, h);

    for (uint dof = 0; dof < elt.data.get_ndof(); dof++) {
        subvector(internal.w1_rhs, SWE::n_dimensions * dof, SWE::n_dimensions) =
            elt.IntegrationPhi(dof, internal.w1_rhs_kernel_at_gp);
    }
}
}
}

#endif
