#ifndef EHDG_GN_PRE_SERIAL_HPP
#define EHDG_GN_PRE_SERIAL_HPP

#include "ehdg_gn_pre_dbath_serial.hpp"

namespace GN {
namespace EHDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  const ProblemInputType& problem_specific_input) {
    SWE_SIM::Problem::preprocessor_serial(discretization, problem_specific_input);

    Problem::initialize_dc_data_serial(discretization.mesh, problem_specific_input);

    uint dc_global_dof_offset = 0;

    Problem::initialize_global_dc_problem_serial(discretization, dc_global_dof_offset);

    discretization.global_data.w1_hat_w1_hat.resize(dc_global_dof_offset, dc_global_dof_offset);
    discretization.global_data.w1_hat_rhs.resize(dc_global_dof_offset);

    Problem::compute_bathymetry_derivatives_serial(discretization);
}
}
}

#endif