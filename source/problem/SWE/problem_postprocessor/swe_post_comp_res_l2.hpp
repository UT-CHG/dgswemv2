#ifndef SWE_POST_COMP_RES_L2_HPP
#define SWE_POST_COMP_RES_L2_HPP

#include "problem/SWE/problem_function_files/swe_true_solution_functions.hpp"

namespace SWE {
template <typename StepperType, typename ElementType>
double compute_residual_L2(const StepperType& stepper, ElementType& elt) {
    double t = stepper.GetTimeAtCurrentStage();

    auto true_q = [t](Point<2>& pt) { return SWE::true_q(t, pt); };

    DynMatrix<double> tmp(SWE::n_variables, columns(elt.data.state[0].q[0]));
    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        row(tmp,var) = elt.data.state[0].q[var];
    }
    return elt.ComputeResidualL2(true_q, tmp);
}
}

#endif