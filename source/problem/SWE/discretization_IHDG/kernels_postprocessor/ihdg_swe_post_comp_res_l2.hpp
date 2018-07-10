#ifndef IHDG_SWE_POST_COMP_RES_L2_HPP
#define IHDG_SWE_POST_COMP_RES_L2_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
template <typename ElementType>
double Problem::compute_residual_L2_kernel(const RKStepper& stepper, ElementType& elt) {
    double t = stepper.GetTimeAtCurrentStage();

    auto true_u = [t](Point<2>& pt) { return SWE::true_u(t, pt); };

    return elt.ComputeResidualL2(true_u, elt.data.state[0].q)[0];
}
}
}

#endif