#ifndef IHDG_GN_PRE_SERIAL_HPP
#define IHDG_GN_PRE_SERIAL_HPP

#include "ihdg_gn_pre_init_data.hpp"
#include "ihdg_gn_pre_init_global_prob.hpp"
#include "ihdg_gn_pre_dbath_serial.hpp"

namespace GN {
namespace IHDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_serial(discretization.mesh, problem_specific_input);

    Problem::initialize_global_problem(discretization);

    Problem::compute_bathymetry_derivatives_serial(discretization);
}
}
}

#endif