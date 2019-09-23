#ifndef EHDG_GN_PROC_SOURCE_HPP
#define EHDG_GN_PROC_SOURCE_HPP

namespace GN {
namespace EHDG {
template <typename ElementType>
void Problem::local_dc_source_kernel(const ESSPRKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();
    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    // at this point h_at_gp, u_at_gp, du_at_gp and ddu_at_gp
    // have been calculated in derivatives kernel

    const auto h  = row(internal.aux_at_gp, SWE::Auxiliaries::h);
    const auto h2 = vec_cw_mult(h, h);
    const auto h3 = vec_cw_mult(h2, h);

    const auto bx = row(internal.dbath_at_gp, GlobalCoord::x);
    const auto by = row(internal.dbath_at_gp, GlobalCoord::y);

    const auto zex = elt.ComputeUgp(row(state.dze, GlobalCoord::x));
    const auto zey = elt.ComputeUgp(row(state.dze, GlobalCoord::y));

    const auto hx = zex + bx;
    const auto hy = zey + by;

    const auto u = row(internal.u_at_gp, GlobalCoord::x);
    const auto v = row(internal.u_at_gp, GlobalCoord::y);

    const auto ux = row(internal.du_at_gp, GN::DU::ux);
    const auto uy = row(internal.du_at_gp, GN::DU::uy);
    const auto vx = row(internal.du_at_gp, GN::DU::vx);
    const auto vy = row(internal.du_at_gp, GN::DU::vy);

    const auto uxx = row(internal.ddu_at_gp, GN::DDU::uxx);
    const auto uxy = row(internal.ddu_at_gp, GN::DDU::uxy);
    const auto uyx = row(internal.ddu_at_gp, GN::DDU::uyx);
    const auto uyy = row(internal.ddu_at_gp, GN::DDU::uyy);
    const auto vxx = row(internal.ddu_at_gp, GN::DDU::vxx);
    const auto vxy = row(internal.ddu_at_gp, GN::DDU::vxy);
    const auto vyx = row(internal.ddu_at_gp, GN::DDU::vyx);
    const auto vyy = row(internal.ddu_at_gp, GN::DDU::vyy);

    const auto bxx = row(internal.ddbath_at_gp, GN::DDBath::bxx);
    const auto bxy = row(internal.ddbath_at_gp, GN::DDBath::bxy);
    const auto byx = row(internal.ddbath_at_gp, GN::DDBath::byx);
    const auto byy = row(internal.ddbath_at_gp, GN::DDBath::byy);

    const auto bxxx = row(internal.dddbath_at_gp, GN::DDDBath::bxxx);
    const auto bxxy = row(internal.dddbath_at_gp, GN::DDDBath::bxxy);
    const auto bxyx = row(internal.dddbath_at_gp, GN::DDDBath::bxyx);
    const auto bxyy = row(internal.dddbath_at_gp, GN::DDDBath::bxyy);
    const auto byxx = row(internal.dddbath_at_gp, GN::DDDBath::byxx);
    const auto byxy = row(internal.dddbath_at_gp, GN::DDDBath::byxy);
    const auto byyx = row(internal.dddbath_at_gp, GN::DDDBath::byyx);
    const auto byyy = row(internal.dddbath_at_gp, GN::DDDBath::byyy);

    const auto c1 = vec_cw_mult(vx, uy) + vec_cw_mult(ux, ux) + vec_cw_mult(ux, vy) + vec_cw_mult(vy, vy);

    const auto c2 = vec_cw_mult(uy, vxx) + vec_cw_mult(vx, uyx) + 2.0 * vec_cw_mult(ux, uxx) + vec_cw_mult(vy, uxx) +
                    vec_cw_mult(ux, vyx) + 2.0 * vec_cw_mult(vy, vyx);

    const auto c3 = vec_cw_mult(uy, vxy) + vec_cw_mult(vx, uyy) + 2.0 * vec_cw_mult(ux, uxy) + vec_cw_mult(vy, uxy) +
                    vec_cw_mult(ux, vyy) + 2.0 * vec_cw_mult(vy, vyy);

    const auto c4 = vec_cw_mult(vec_cw_mult(u, u), bxx) + vec_cw_mult(vec_cw_mult(u, v), bxy) +
                    vec_cw_mult(vec_cw_mult(u, v), byx) + vec_cw_mult(vec_cw_mult(v, v), byy);

    const auto c5 = 2.0 * vec_cw_mult(vec_cw_mult(u, ux), bxx) + vec_cw_mult(vec_cw_mult(u, u), bxxx) +
                    vec_cw_mult(vec_cw_mult(ux, v), bxy) + vec_cw_mult(vec_cw_mult(u, vx), bxy) +
                    vec_cw_mult(vec_cw_mult(u, v), bxyx) + vec_cw_mult(vec_cw_mult(ux, v), byx) +
                    vec_cw_mult(vec_cw_mult(u, vx), byx) + vec_cw_mult(vec_cw_mult(u, v), byxx) +
                    2.0 * vec_cw_mult(vec_cw_mult(v, vx), byy) + vec_cw_mult(vec_cw_mult(v, v), byyx);

    const auto c6 = 2.0 * vec_cw_mult(vec_cw_mult(u, uy), bxx) + vec_cw_mult(vec_cw_mult(u, u), bxxy) +
                    vec_cw_mult(vec_cw_mult(uy, v), bxy) + vec_cw_mult(vec_cw_mult(u, vy), bxy) +
                    vec_cw_mult(vec_cw_mult(u, v), bxyy) + vec_cw_mult(vec_cw_mult(uy, v), byx) +
                    vec_cw_mult(vec_cw_mult(u, vy), byx) + vec_cw_mult(vec_cw_mult(u, v), byxy) +
                    2.0 * vec_cw_mult(vec_cw_mult(v, vy), byy) + vec_cw_mult(vec_cw_mult(v, v), byyy);

    row(internal.w1_rhs_kernel_at_gp, GlobalCoord::x) =
        2.0 * vec_cw_mult(vec_cw_mult(h2, c1), hx) + 2.0 / 3.0 * vec_cw_mult(h3, c2) +
        vec_cw_mult(vec_cw_mult(h, c4), hx) + 1.0 / 2.0 * vec_cw_mult(h2, c5) + vec_cw_mult(vec_cw_mult(h2, c1), bx) +
        vec_cw_mult(vec_cw_mult(h, c4), bx) + Global::g / NDParameters::alpha * vec_cw_mult(zex, h);

    row(internal.w1_rhs_kernel_at_gp, GlobalCoord::y) =
        2.0 * vec_cw_mult(vec_cw_mult(h2, c1), hy) + 2.0 / 3.0 * vec_cw_mult(h3, c3) +
        vec_cw_mult(vec_cw_mult(h, c4), hy) + 1.0 / 2.0 * vec_cw_mult(h2, c6) + vec_cw_mult(vec_cw_mult(h2, c1), by) +
        vec_cw_mult(vec_cw_mult(h, c4), by) + Global::g / NDParameters::alpha * vec_cw_mult(zey, h);

    for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
        subvector(internal.w1_rhs, GN::n_dimensions * dof, GN::n_dimensions) =
            elt.IntegrationPhi(dof, internal.w1_rhs_kernel_at_gp);
    }
}
}
}

#endif
