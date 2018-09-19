#ifndef IHDG_SWE_PRE_SERIAL_HPP
#define IHDG_SWE_PRE_SERIAL_HPP

#include "ihdg_swe_pre_init_data.hpp"

namespace SWE {
namespace IHDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_serial(discretization.mesh, problem_specific_input);

    uint global_dof_offset = 0;

    Problem::initialize_global_problem_serial(discretization, global_dof_offset);
}
}
}

#endif