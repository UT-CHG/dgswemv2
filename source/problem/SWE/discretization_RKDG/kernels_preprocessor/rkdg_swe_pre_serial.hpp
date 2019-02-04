#ifndef RKDG_SWE_PRE_SERIAL_HPP
#define RKDG_SWE_PRE_SERIAL_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace RKDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  ProblemGlobalDataType& global_data,
                                  const ProblemInputType& problem_specific_input) {
    initialize_data_serial(discretization.mesh, problem_specific_input);
}
}
}

#endif