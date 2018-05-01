#ifndef SWE_POST_COMP_RES_L2_HPP
#define SWE_POST_COMP_RES_L2_HPP

#include "../../../general_definitions.hpp"

namespace SWE {
template <typename ElementType>
double Problem::compute_residual_L2_kernel(const Stepper& stepper, ElementType& elt) {
    double t = stepper.GetTimeAtCurrentStage();

    auto true_ze = [t](Point<2>& pt) { return SWE::true_ze(t, pt); };

    return elt.ComputeResidualL2(true_ze, elt.data.state[0].ze);
}
}

#endif