#ifndef IHDG_SWE_PROC_UPDATE_HPP
#define IHDG_SWE_PROC_UPDATE_HPP

namespace SWE {
namespace IHDG {
template <typename ElementType>
void Problem::update_kernel(const RKStepper& stepper, ElementType& elt) {}

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

    for (auto& u_mode : row(state.q, SWE::Variables::ze)) {
        if (std::isnan(u_mode)) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& u_mode : row(state.q, SWE::Variables::qx)) {
        if (std::isnan(u_mode)) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (auto& u_mode : row(state.q, SWE::Variables::qy)) {
        if (std::isnan(u_mode)) {
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
