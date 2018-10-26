#ifndef SWE_POST_SCRUTINIZE_HPP
#define SWE_POST_SCRUTINIZE_HPP

namespace SWE {
template <typename StepperType, typename ElementType>
bool scrutinize_solution(const StepperType& stepper, ElementType& elt) {
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
        if (std::isnan(state.q(SWE::Variables::qy, dof))) {
            std::cerr << "Error: found isnan qy at Element " << elt.GetID();
            std::cerr << "       At stage: " << stage << "\n";

            return true;
        }
    }

    return false;
}
}

#endif
