#ifndef EHDG_SWE_PRE_SERIAL_HPP
#define EHDG_SWE_PRE_SERIAL_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace EHDG {
template <typename ProblemType>
void Problem::preprocessor_serial(HDGDiscretization<ProblemType>& discretization,
                                  typename ProblemType::ProblemGlobalDataType& global_data,
                                  const typename ProblemType::ProblemStepperType& stepper,
                                  const typename ProblemType::ProblemInputType& problem_specific_input) {
    initialize_data_serial(discretization.mesh, problem_specific_input);

    Problem::initialize_global_problem_serial(discretization);
}
}
}

#endif