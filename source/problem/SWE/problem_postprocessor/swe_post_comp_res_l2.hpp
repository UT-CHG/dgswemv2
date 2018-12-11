#ifndef SWE_POST_COMP_RES_L2_HPP
#define SWE_POST_COMP_RES_L2_HPP

#include "problem/SWE/problem_function_files/fwd.hpp"

namespace SWE {
template <typename StepperType, typename ElementType>
double compute_residual_L2(const StepperType& stepper, ElementType& elt) {
    double t = stepper.GetTimeAtCurrentStage();

    auto true_q = [t](Point<2>& pt) { return SWE::true_q(t, pt); };

    return elt.ComputeResidualL2(true_q, elt.data.state[0].q);
}
}

#endif