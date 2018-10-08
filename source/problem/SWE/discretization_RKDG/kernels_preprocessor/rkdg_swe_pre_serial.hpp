#ifndef RKDG_SWE_PRE_SERIAL_HPP
#define RKDG_SWE_PRE_SERIAL_HPP

namespace SWE {
namespace RKDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_serial(discretization.mesh, problem_specific_input);
}
}
}

#endif