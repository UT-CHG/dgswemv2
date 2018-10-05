#ifndef EHDG_GN_PRE_INIT_PROB_PARAMS_HPP
#define EHDG_GN_PRE_INIT_PROB_PARAMS_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_problem_parameters(const ProblemInputType& problem_specific_input) {
    GN::Global::g         = problem_specific_input.g;
    GN::Global::rho_air   = problem_specific_input.rho_air;
    GN::Global::rho_water = problem_specific_input.rho_water;

    // specify forcing terms
    if (problem_specific_input.function_source.type != SWE::FunctionSourceType::None) {
        SWE::SourceTerms::function_source = true;
    }

    if (problem_specific_input.bottom_friction.type != SWE::BottomFrictionType::None) {
        SWE::SourceTerms::bottom_friction = true;
        SWE::SourceTerms::Cf              = problem_specific_input.bottom_friction.coefficient;
    }

    if (problem_specific_input.meteo_forcing.type != SWE::MeteoForcingType::None) {
        SWE::SourceTerms::meteo_forcing = true;
    }

    if (problem_specific_input.tide_potential.type != SWE::TidePotentialType::None) {
        SWE::SourceTerms::tide_potential = true;
    }

    if (problem_specific_input.coriolis.type != SWE::CoriolisType::None) {
        SWE::SourceTerms::coriolis = true;
    }
}
}
}

#endif