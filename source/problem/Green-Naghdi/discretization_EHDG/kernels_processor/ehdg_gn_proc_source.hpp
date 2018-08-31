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
void Problem::local_dc_source_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    // at this point h_at_gp, u_at_gp, du_at_gp and ddu_at_gp
    // have been calculated in derivatives kernel

    auto h  = row(internal.aux_at_gp, GN::Auxiliaries::h);
    auto h2 = cwise_multiplication(h, h);
    auto h3 = cwise_multiplication(h2, h);

    auto bx = row(internal.dbath_at_gp, GlobalCoord::x);
    auto by = row(internal.dbath_at_gp, GlobalCoord::y);

    auto zex = elt.ComputeUgp(row(state.dze, GlobalCoord::x));
    auto zey = elt.ComputeUgp(row(state.dze, GlobalCoord::y));

    auto hx = zex + bx;
    auto hy = zey + by;

    auto u = row(internal.u_at_gp, GlobalCoord::x);
    auto v = row(internal.u_at_gp, GlobalCoord::y);

    auto ux = row(internal.du_at_gp, GN::DU::ux);
    auto uy = row(internal.du_at_gp, GN::DU::uy);
    auto vx = row(internal.du_at_gp, GN::DU::vx);
    auto vy = row(internal.du_at_gp, GN::DU::vy);

    auto uxx = row(internal.ddu_at_gp, GN::DDU::uxx);
    auto uxy = row(internal.ddu_at_gp, GN::DDU::uxy);
    auto uyx = row(internal.ddu_at_gp, GN::DDU::uyx);
    auto uyy = row(internal.ddu_at_gp, GN::DDU::uyy);
    auto vxx = row(internal.ddu_at_gp, GN::DDU::vxx);
    auto vxy = row(internal.ddu_at_gp, GN::DDU::vxy);
    auto vyx = row(internal.ddu_at_gp, GN::DDU::vyx);
    auto vyy = row(internal.ddu_at_gp, GN::DDU::vyy);

    auto bxx = row(internal.ddbath_at_gp, GN::DDBath::bxx);
    auto bxy = row(internal.ddbath_at_gp, GN::DDBath::bxy);
    auto byx = row(internal.ddbath_at_gp, GN::DDBath::byx);
    auto byy = row(internal.ddbath_at_gp, GN::DDBath::byy);

    auto bxxx = row(internal.dddbath_at_gp, GN::DDDBath::bxxx);
    auto bxxy = row(internal.dddbath_at_gp, GN::DDDBath::bxxy);
    auto bxyx = row(internal.dddbath_at_gp, GN::DDDBath::bxyx);
    auto bxyy = row(internal.dddbath_at_gp, GN::DDDBath::bxyy);
    auto byxx = row(internal.dddbath_at_gp, GN::DDDBath::byxx);
    auto byxy = row(internal.dddbath_at_gp, GN::DDDBath::byxy);
    auto byyx = row(internal.dddbath_at_gp, GN::DDDBath::byyx);
    auto byyy = row(internal.dddbath_at_gp, GN::DDDBath::byyy);

    auto c1 = cwise_multiplication(vx, uy) + cwise_multiplication(ux, ux) + cwise_multiplication(ux, vy) +
              cwise_multiplication(vy, vy);

    auto c2 = cwise_multiplication(uy, vxx) + cwise_multiplication(vx, uyx) + 2.0 * cwise_multiplication(ux, uxx) +
              cwise_multiplication(vy, uxx) + cwise_multiplication(ux, vyx) + 2.0 * cwise_multiplication(vy, vyx);

    auto c3 = cwise_multiplication(uy, vxy) + cwise_multiplication(vx, uyy) + 2.0 * cwise_multiplication(ux, uxy) +
              cwise_multiplication(vy, uxy) + cwise_multiplication(ux, vyy) + 2.0 * cwise_multiplication(vy, vyy);

    auto c4 =
        cwise_multiplication(cwise_multiplication(u, u), bxx) + cwise_multiplication(cwise_multiplication(u, v), bxy) +
        cwise_multiplication(cwise_multiplication(u, v), byx) + cwise_multiplication(cwise_multiplication(v, v), byy);

    auto c5 = 2.0 * cwise_multiplication(cwise_multiplication(u, ux), bxx) +
              cwise_multiplication(cwise_multiplication(u, u), bxxx) +
              cwise_multiplication(cwise_multiplication(ux, v), bxy) +
              cwise_multiplication(cwise_multiplication(u, vx), bxy) +
              cwise_multiplication(cwise_multiplication(u, v), bxyx) +
              cwise_multiplication(cwise_multiplication(ux, v), byx) +
              cwise_multiplication(cwise_multiplication(u, vx), byx) +
              cwise_multiplication(cwise_multiplication(u, v), byxx) +
              2.0 * cwise_multiplication(cwise_multiplication(v, vx), byy) +
              cwise_multiplication(cwise_multiplication(v, v), byyx);

    auto c6 = 2.0 * cwise_multiplication(cwise_multiplication(u, uy), bxx) +
              cwise_multiplication(cwise_multiplication(u, u), bxxy) +
              cwise_multiplication(cwise_multiplication(uy, v), bxy) +
              cwise_multiplication(cwise_multiplication(u, vy), bxy) +
              cwise_multiplication(cwise_multiplication(u, v), bxyy) +
              cwise_multiplication(cwise_multiplication(uy, v), byx) +
              cwise_multiplication(cwise_multiplication(u, vy), byx) +
              cwise_multiplication(cwise_multiplication(u, v), byxy) +
              2.0 * cwise_multiplication(cwise_multiplication(v, vy), byy) +
              cwise_multiplication(cwise_multiplication(v, v), byyy);

    row(internal.w1_rhs_kernel_at_gp, GlobalCoord::x) =
        2.0 * cwise_multiplication(cwise_multiplication(h2, c1), hx) + 2.0 / 3.0 * cwise_multiplication(h3, c2) +
        cwise_multiplication(cwise_multiplication(h, c4), hx) + 1.0 / 2.0 * cwise_multiplication(h2, c5) +
        cwise_multiplication(cwise_multiplication(h2, c1), bx) + cwise_multiplication(cwise_multiplication(h, c4), bx) +
        Global::g / NDParameters::alpha * cwise_multiplication(zex, h);

    row(internal.w1_rhs_kernel_at_gp, GlobalCoord::y) =
        2.0 * cwise_multiplication(cwise_multiplication(h2, c1), hy) + 2.0 / 3.0 * cwise_multiplication(h3, c3) +
        cwise_multiplication(cwise_multiplication(h, c4), hy) + 1.0 / 2.0 * cwise_multiplication(h2, c6) +
        cwise_multiplication(cwise_multiplication(h2, c1), by) + cwise_multiplication(cwise_multiplication(h, c4), by) +
        Global::g / NDParameters::alpha * cwise_multiplication(zey, h);

    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        subvector(internal.w1_rhs, GN::n_dimensions * dof, GN::n_dimensions) =
            elt.IntegrationPhi(dof, internal.w1_rhs_kernel_at_gp);
    }
}
}
}

#endif
