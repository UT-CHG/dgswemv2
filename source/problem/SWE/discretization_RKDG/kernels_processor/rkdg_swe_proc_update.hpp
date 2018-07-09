#ifndef RKDG_SWE_PROC_UPDATE_HPP
#define RKDG_SWE_PROC_UPDATE_HPP

namespace SWE {
namespace RKDG {
template <typename ElementType>
void Problem::update_kernel(const RKStepper& stepper, ElementType& elt) {
    const uint stage = stepper.GetStage();
    const double dt  = stepper.GetDT();

    auto& state      = elt.data.state;
    auto& curr_state = elt.data.state[stage];
    auto& next_state = elt.data.state[stage + 1];

    elt.ApplyMinv(curr_state.rhs, curr_state.solution);

    std::fill(next_state.q.begin(), next_state.q.end(), 0.0);

    for (uint s = 0; s <= stage; ++s) {
        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            next_state.q[dof] +=
                stepper.ark[stage][s] * state[s].q[dof] + dt * stepper.brk[stage][s] * state[s].solution[dof];
        }
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

    for (auto& u_mode : state.q) {
        if (std::isnan(u_mode[SWE::Variables::ze])) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& u_mode : state.q) {
        if (std::isnan(u_mode[SWE::Variables::qx])) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& u_mode : state.q) {
        if (std::isnan(u_mode[SWE::Variables::qy])) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& rhs_mode : state.rhs) {
        if (std::isnan(rhs_mode[SWE::Variables::ze])) {
            std::cerr << "Error: found isnan rhs_ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& rhs_mode : state.rhs) {
        if (std::isnan(rhs_mode[SWE::Variables::qx])) {
            std::cerr << "Error: found isnan rhs_qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& rhs_mode : state.rhs) {
        if (std::isnan(rhs_mode[SWE::Variables::qy])) {
            std::cerr << "Error: found isnan rhs_qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    return false;
}
}
}

#endif
