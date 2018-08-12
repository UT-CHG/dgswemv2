#ifndef EHDG_GN_PRE_SERIAL_HPP
#define EHDG_GN_PRE_SERIAL_HPP

#include "ehdg_gn_pre_init_data.hpp"

namespace GN {
namespace EHDG {
void Problem::serial_preprocessor_kernel(ProblemDiscretizationType& discretization,
                                         const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_kernel(discretization.mesh, problem_specific_input);

    Problem::initialize_global_problem(discretization);
}
}
}

#endif