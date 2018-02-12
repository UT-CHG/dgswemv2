#ifndef SWE_PRE_INIT_PROB_PARAMS_HPP
#define SWE_PRE_INIT_PROB_PARAMS_HPP

namespace SWE {
void Problem::initialize_problem_parameters(const ProblemInputType& problem_specific_input) {
    // specify forcing terms
    SWE::Global::g = problem_specific_input.g;
    SWE::Global::h_o = problem_specific_input.h_o;

    if (problem_specific_input.function_source.type != SWE::FunctionSourceType::None) {
        SWE::SourceTerms::function_source = true;
    }

    if (problem_specific_input.bottom_friction.type != SWE::BottomFrictionType::None) {
        SWE::SourceTerms::bottom_friction = true;
        SWE::Global::Cf = problem_specific_input.bottom_friction.coefficient;
    }

    if (problem_specific_input.meteo_forcing.type != SWE::MeteoForcingType::None) {
        SWE::SourceTerms::meteo_forcing = true;
    }

    if (problem_specific_input.tidal_potential.type != SWE::TidalPotentialType::None) {
        SWE::SourceTerms::tidal_potential = false;
    }

    if (problem_specific_input.coriolis.type != SWE::CoriolisType::None) {
        SWE::SourceTerms::coriolis = false;
    }
}
}

#endif