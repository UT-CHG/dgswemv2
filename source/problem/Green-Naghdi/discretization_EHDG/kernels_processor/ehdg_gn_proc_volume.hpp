#ifndef EHDG_GN_PROC_VOLUME_HPP
#define EHDG_GN_PROC_VOLUME_HPP

namespace GN {
namespace EHDG {
template <typename ElementType>
void Problem::local_swe_volume_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    row(internal.aux_at_gp, GN::Auxiliaries::h) =
        row(internal.q_at_gp, GN::Variables::ze) + row(internal.aux_at_gp, GN::Auxiliaries::bath);

    auto u = cwise_division(row(internal.q_at_gp, GN::Variables::qx), row(internal.aux_at_gp, GN::Auxiliaries::h));
    auto v = cwise_division(row(internal.q_at_gp, GN::Variables::qy), row(internal.aux_at_gp, GN::Auxiliaries::h));

    auto uuh = cwise_multiplication(u, row(internal.q_at_gp, GN::Variables::qx));
    auto vvh = cwise_multiplication(v, row(internal.q_at_gp, GN::Variables::qy));
    auto uvh = cwise_multiplication(u, row(internal.q_at_gp, GN::Variables::qy));
    auto pe  = Global::g * (0.5 * cwise_multiplication(row(internal.q_at_gp, GN::Variables::ze),
                                                      row(internal.q_at_gp, GN::Variables::ze)) +
                           cwise_multiplication(row(internal.q_at_gp, GN::Variables::ze),
                                                row(internal.aux_at_gp, GN::Auxiliaries::bath)));

    row(internal.Fx_at_gp, GN::Variables::ze) = row(internal.q_at_gp, GN::Variables::qx);
    row(internal.Fx_at_gp, GN::Variables::qx) = uuh + pe;
    row(internal.Fx_at_gp, GN::Variables::qy) = uvh;

    row(internal.Fy_at_gp, GN::Variables::ze) = row(internal.q_at_gp, GN::Variables::qy);
    row(internal.Fy_at_gp, GN::Variables::qx) = uvh;
    row(internal.Fy_at_gp, GN::Variables::qy) = vvh + pe;

    state.rhs =
        elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) + elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
}
}
}

#endif
