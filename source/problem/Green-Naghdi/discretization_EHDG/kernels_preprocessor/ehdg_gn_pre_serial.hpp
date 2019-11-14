#ifndef EHDG_GN_PRE_SERIAL_HPP
#define EHDG_GN_PRE_SERIAL_HPP

#include "problem/Green-Naghdi/problem_preprocessor/gn_pre_init_data.hpp"
#include "ehdg_gn_pre_dbath_serial.hpp"

namespace GN {
namespace EHDG {
void Problem::preprocessor_serial(ProblemDiscretizationType& discretization,
                                  ProblemGlobalDataType& global_data,
                                  const ProblemStepperType& stepper,
                                  const ProblemInputType& problem_specific_input) {
    SWE_SIM::Problem::preprocessor_serial(
        discretization, global_data, stepper.GetFirstStepper(), problem_specific_input);

    GN::initialize_data_serial(discretization.mesh);

    uint dc_global_dof_offset = 0;
    Problem::initialize_global_dc_problem_serial(discretization, dc_global_dof_offset);
    global_data.w1_hat_w1_hat.resize(dc_global_dof_offset, dc_global_dof_offset);
    global_data.w1_hat_rhs.resize(dc_global_dof_offset);

    Problem::compute_bathymetry_derivatives_serial(discretization, global_data);

    uint n_stages = stepper.GetFirstStepper().GetNumStages() > stepper.GetSecondStepper().GetNumStages()
                        ? stepper.GetFirstStepper().GetNumStages()
                        : stepper.GetSecondStepper().GetNumStages();

    discretization.mesh.CallForEachElement([n_stages](auto& elt) { elt.data.resize(n_stages + 1); });
}
}
}

#endif