#ifndef IHDG_SWE_PROC_UPDATE_HPP
#define IHDG_SWE_PROC_UPDATE_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename ElementType>
void Problem::update_kernel(const StepperType& stepper, ElementType& elt) {}

template <typename StepperType, typename ElementType>
void Problem::swap_states_kernel(const StepperType& stepper, ElementType& elt) {
    uint n_stages = stepper.GetNumStages();
    auto& state   = elt.data.state;

    std::swap(state[0].q, state[n_stages].q);
}

template <typename StepperType, typename ElementType>
bool Problem::scrutinize_solution_kernel(const StepperType& stepper, ElementType& elt) {
    uint stage = stepper.GetStage();

    auto& state = elt.data.state[stage + 1];

    uint ndof = elt.data.get_ndof();

    for (uint dof = 0; dof < ndof; ++dof) {
        if (std::isnan(state.q(SWE::Variables::ze, dof))) {
            std::cerr << "Error: found isnan ze at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        if (std::isnan(state.q(SWE::Variables::qx, dof))) {
            std::cerr << "Error: found isnan qx at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        if (std::isnan(state.q(SWE::Variables::qx, dof))) {
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
