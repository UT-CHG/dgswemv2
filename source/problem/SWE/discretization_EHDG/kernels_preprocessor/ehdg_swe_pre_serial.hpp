#ifndef EHDG_SWE_PRE_SERIAL_HPP
#define EHDG_SWE_PRE_SERIAL_HPP

#include "ehdg_swe_pre_init_data.hpp"

namespace SWE {
namespace EHDG {
void Problem::serial_preprocessor_kernel(ProblemDiscretizationType& discretization,
                                         const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_serial_kernel(discretization.mesh, problem_specific_input);

    Problem::initialize_global_problem(discretization);
}
}
}

#endif