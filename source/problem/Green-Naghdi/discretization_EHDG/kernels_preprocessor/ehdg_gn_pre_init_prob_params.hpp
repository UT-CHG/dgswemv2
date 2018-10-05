#ifndef EHDG_GN_PRE_INIT_PROB_PARAMS_HPP
#define EHDG_GN_PRE_INIT_PROB_PARAMS_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_problem_parameters(const ProblemInputType& problem_specific_input) {
    SWE::EHDG::Problem::initialize_problem_parameters(problem_specific_input);

    GN::Global::g         = problem_specific_input.g;
    GN::Global::rho_air   = problem_specific_input.rho_air;
    GN::Global::rho_water = problem_specific_input.rho_water;
}
}
}

#endif