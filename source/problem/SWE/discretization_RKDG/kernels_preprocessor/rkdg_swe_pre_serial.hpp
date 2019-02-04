#ifndef RKDG_SWE_PRE_SERIAL_HPP
#define RKDG_SWE_PRE_SERIAL_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace RKDG {
template <typename SerialSimType>
void Problem::preprocessor_serial(SerialSimType* sim) {
    auto& discretization         = sim->discretization;
    auto& problem_specific_input = sim->problem_input;

    initialize_data_serial(discretization.mesh, problem_specific_input);
}
}
}

#endif