#ifndef IHDG_SWE_PRE_SERIAL_HPP
#define IHDG_SWE_PRE_SERIAL_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace IHDG {
template <typename ProblemType>
void Problem::preprocessor_serial(HDGDiscretization<ProblemType>& discretization,
                                  const ProblemInputType& problem_specific_input) {
    initialize_data_serial(discretization.mesh, problem_specific_input);

    uint global_dof_offset = 0;

    Problem::initialize_global_problem_serial(discretization, global_dof_offset);

    discretization.global_data.delta_hat_global.resize(global_dof_offset, global_dof_offset);
    discretization.global_data.rhs_global.resize(global_dof_offset);
}
}
}

#endif