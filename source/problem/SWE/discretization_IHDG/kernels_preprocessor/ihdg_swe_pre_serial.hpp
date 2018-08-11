#ifndef IHDG_SWE_PRE_SERIAL_HPP
#define IHDG_SWE_PRE_SERIAL_HPP

#include "ihdg_swe_pre_init_data.hpp"

namespace SWE {
namespace IHDG {
void Problem::serial_preprocessor_kernel(ProblemDiscretizationType& discretization,
                                         const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_kernel(discretization.mesh, problem_specific_input);

    Problem::initialize_global_problem(discretization);
}
}
}

#endif