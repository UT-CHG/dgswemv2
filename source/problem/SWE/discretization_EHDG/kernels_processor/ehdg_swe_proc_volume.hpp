#ifndef EHDG_SWE_PROC_VOLUME_HPP
#define EHDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace EHDG {
template <typename StepperType, typename ElementType>
void Problem::local_volume_kernel(const StepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    internal.q_at_gp = elt.ComputeUgp(state.q);

    row(internal.aux_at_gp, SWE::Auxiliaries::h) =
        row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

    auto u = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(internal.q_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.q_at_gp, SWE::Variables::ze)) +
               vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.aux_at_gp, SWE::Auxiliaries::bath)));

    row(internal.Fx_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qx);
    row(internal.Fx_at_gp, SWE::Variables::qx) = uuh + pe;
    row(internal.Fx_at_gp, SWE::Variables::qy) = uvh;

    row(internal.Fy_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qy);
    row(internal.Fy_at_gp, SWE::Variables::qx) = uvh;
    row(internal.Fy_at_gp, SWE::Variables::qy) = vvh + pe;

    state.rhs =
        elt.IntegrationDPhi(GlobalCoord::x, internal.Fx_at_gp) + elt.IntegrationDPhi(GlobalCoord::y, internal.Fy_at_gp);
}
}
}

#endif
