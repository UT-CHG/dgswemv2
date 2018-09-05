#ifndef IHDG_GN_POST_COMP_RES_L2_HPP
#define IHDG_GN_POST_COMP_RES_L2_HPP

#include "general_definitions.hpp"

namespace GN {
namespace IHDG {
template <typename ElementType>
double Problem::compute_residual_L2(const RKStepper& stepper, ElementType& elt) {
    double t = stepper.GetTimeAtCurrentStage();

    auto true_u = [t](Point<2>& pt) { return GN::true_u(t, pt); };

    return elt.ComputeResidualL2(true_u, elt.data.state[0].q);
}
}
}

#endif