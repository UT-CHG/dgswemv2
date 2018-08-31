#ifndef IHDG_GN_PROC_UPDATE_HPP
#define IHDG_GN_PROC_UPDATE_HPP

namespace GN {
namespace IHDG {
template <typename ElementType>
void Problem::update_kernel(const RKStepper& stepper, ElementType& elt) {
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
void Problem::swap_states_kernel(const RKStepper& stepper, ElementType& elt) {
    uint n_stages = stepper.GetNumStages();
    auto& state   = elt.data.state;

    std::swap(state[0].q, state[n_stages].q);
}

template <typename ElementType>
bool Problem::scrutinize_solution_kernel(const RKStepper& stepper, ElementType& elt) {
    uint stage = stepper.GetStage();

    auto& state = elt.data.state[stage + 1];

    uint ndof = elt.data.get_ndof();

    for (uint dof = 0; dof < ndof; ++dof) {
        if (std::isnan(state.q(GN::Variables::ze, dof))) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        if (std::isnan(state.q(GN::Variables::qx, dof))) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        if (std::isnan(state.q(GN::Variables::qx, dof))) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    return false;
}
}
}

#endif