#ifndef RKDG_SWE_PROC_VOLUME_HPP
#define RKDG_SWE_PROC_VOLUME_HPP

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::volume_kernel(const ProblemStepperType& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& wd_state = elt.data.wet_dry_state;
    auto& state    = elt.data.state[stage];

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        set_constant(state.rhs[var], 0.0);
    }

    if (wd_state.wet) {
        auto& internal = elt.data.internal;

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            row(internal.q_at_gp,var) = elt.ComputeUgp(state.q[var]);
        }

        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        auto u = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
        auto v = vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

        auto uuh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qx));
        auto vvh = vec_cw_mult(v, row(internal.q_at_gp, SWE::Variables::qy));
        auto uvh = vec_cw_mult(u, row(internal.q_at_gp, SWE::Variables::qy));
        auto pe =
            Global::g *
            (0.5 * vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.q_at_gp, SWE::Variables::ze)) +
             vec_cw_mult(row(internal.q_at_gp, SWE::Variables::ze), row(internal.aux_at_gp, SWE::Auxiliaries::bath)));

        row(internal.Fx_at_gp, SWE::Variables::ze) =
            vec_cw_mult(row(internal.aux_at_gp, SWE::Auxiliaries::sp), row(internal.q_at_gp, SWE::Variables::qx));
        row(internal.Fx_at_gp, SWE::Variables::qx) =
            vec_cw_mult(row(internal.aux_at_gp, SWE::Auxiliaries::sp), uuh + pe);
        row(internal.Fx_at_gp, SWE::Variables::qy) = vec_cw_mult(row(internal.aux_at_gp, SWE::Auxiliaries::sp), uvh);

        row(internal.Fy_at_gp, SWE::Variables::ze) = row(internal.q_at_gp, SWE::Variables::qy);
        row(internal.Fy_at_gp, SWE::Variables::qx) = uvh;
        row(internal.Fy_at_gp, SWE::Variables::qy) = vvh + pe;

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            state.rhs[var] = elt.IntegrationDPhi(GlobalCoord::x, row(internal.Fx_at_gp,var)) +
                elt.IntegrationDPhi(GlobalCoord::y, row(internal.Fy_at_gp,var));
        }
    }
}
}
}

#endif
