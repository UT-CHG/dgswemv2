#ifndef EHDG_GN_PROC_UPDATE_HPP
#define EHDG_GN_PROC_UPDATE_HPP

namespace GN {
namespace EHDG {
template <typename ElementType>
void Problem::dc_update_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();
    const double dt  = stepper.GetDT();

    auto& state      = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    curr_state.solution = elt.ApplyMinv(curr_state.rhs);

    set_constant(next_state.q, 0.0);

    for (uint s = 0; s <= stage; ++s) {
        next_state.q += stepper.ark[stage][s] * state[s].q + dt * stepper.brk[stage][s] * state[s].solution;
    }
}

template <typename ElementType>
void Problem::dispersive_correction_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();

    auto& state    = elt.data.state[stage];
    auto& internal = elt.data.internal;

    auto h = row(internal.aux_at_gp, SWE::Auxiliaries::h);

    auto dze_dx = elt.ComputeUgp(row(state.dze, GlobalCoord::x));
    auto dze_dy = elt.ComputeUgp(row(state.dze, GlobalCoord::y));

    row(internal.source_at_gp, SWE::Variables::qx) = Global::g / NDParameters::alpha * vec_cw_mult(dze_dx, h);
    row(internal.source_at_gp, SWE::Variables::qy) = Global::g / NDParameters::alpha * vec_cw_mult(dze_dy, h);

    row(internal.source_at_gp, SWE::Variables::qx) -= elt.ComputeUgp(row(state.w1, GlobalCoord::x));
    row(internal.source_at_gp, SWE::Variables::qy) -= elt.ComputeUgp(row(state.w1, GlobalCoord::y));

    row(state.rhs, SWE::Variables::qx) = elt.IntegrationPhi(row(internal.source_at_gp, SWE::Variables::qx));
    row(state.rhs, SWE::Variables::qy) = elt.IntegrationPhi(row(internal.source_at_gp, SWE::Variables::qy));

    row(state.solution, SWE::Variables::qx) = elt.ApplyMinv(row(state.rhs, SWE::Variables::qx));
    row(state.solution, SWE::Variables::qy) = elt.ApplyMinv(row(state.rhs, SWE::Variables::qy));

    row(state.q, SWE::Variables::qx) += stepper.GetDT() * row(state.solution, SWE::Variables::qx);
    row(state.q, SWE::Variables::qy) += stepper.GetDT() * row(state.solution, SWE::Variables::qy);
}

template <typename ElementType>
void Problem::swap_states_kernel(const RKStepper& stepper, ElementType& elt) {
    uint n_stages = stepper.GetNumStages();
    auto& state   = elt.data.state;

    std::swap(state[0].q, state[n_stages].q);
}
}
}

#endif
