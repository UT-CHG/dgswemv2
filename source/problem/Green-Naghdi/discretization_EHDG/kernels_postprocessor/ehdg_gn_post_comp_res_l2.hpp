#ifndef EHDG_GN_POST_COMP_RES_L2_HPP
#define EHDG_GN_POST_COMP_RES_L2_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
template <typename ElementType>
double Problem::compute_residual_L2(const RKStepper& stepper, ElementType& elt) {
    return SWE::EHDG::Problem::compute_residual_L2(stepper, elt);
}
}
}

#endif