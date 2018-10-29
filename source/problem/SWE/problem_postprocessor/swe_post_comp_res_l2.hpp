#ifndef SWE_POST_COMP_RES_L2_HPP
#define SWE_POST_COMP_RES_L2_HPP

#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

namespace SWE {
template <typename StepperType, typename ElementType>
double compute_residual_L2(const StepperType& stepper, ElementType& elt) {
    double t = stepper.GetTimeAtCurrentStage();

    auto true_u = [t](Point<2>& pt) { return SWE::true_u(t, pt); };

    return elt.ComputeResidualL2(true_u, elt.data.state[0].q);
}
}

#endif