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

    auto h  = row(internal.aux_at_gp, SWE::Auxiliaries::h);
    auto h2 = vec_cw_mult(h, h);
    auto h3 = vec_cw_mult(h2, h);

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

    auto c1 = vec_cw_mult(vx, uy) + vec_cw_mult(ux, ux) + vec_cw_mult(ux, vy) + vec_cw_mult(vy, vy);

    auto c2 = vec_cw_mult(uy, vxx) + vec_cw_mult(vx, uyx) + 2.0 * vec_cw_mult(ux, uxx) + vec_cw_mult(vy, uxx) +
              vec_cw_mult(ux, vyx) + 2.0 * vec_cw_mult(vy, vyx);

    auto c3 = vec_cw_mult(uy, vxy) + vec_cw_mult(vx, uyy) + 2.0 * vec_cw_mult(ux, uxy) + vec_cw_mult(vy, uxy) +
              vec_cw_mult(ux, vyy) + 2.0 * vec_cw_mult(vy, vyy);

    auto c4 = vec_cw_mult(vec_cw_mult(u, u), bxx) + vec_cw_mult(vec_cw_mult(u, v), bxy) +
              vec_cw_mult(vec_cw_mult(u, v), byx) + vec_cw_mult(vec_cw_mult(v, v), byy);

    auto c5 = 2.0 * vec_cw_mult(vec_cw_mult(u, ux), bxx) + vec_cw_mult(vec_cw_mult(u, u), bxxx) +
              vec_cw_mult(vec_cw_mult(ux, v), bxy) + vec_cw_mult(vec_cw_mult(u, vx), bxy) +
              vec_cw_mult(vec_cw_mult(u, v), bxyx) + vec_cw_mult(vec_cw_mult(ux, v), byx) +
              vec_cw_mult(vec_cw_mult(u, vx), byx) + vec_cw_mult(vec_cw_mult(u, v), byxx) +
              2.0 * vec_cw_mult(vec_cw_mult(v, vx), byy) + vec_cw_mult(vec_cw_mult(v, v), byyx);

    auto c6 = 2.0 * vec_cw_mult(vec_cw_mult(u, uy), bxx) + vec_cw_mult(vec_cw_mult(u, u), bxxy) +
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
