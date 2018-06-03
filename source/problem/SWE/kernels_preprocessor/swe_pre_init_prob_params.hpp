#ifndef SWE_PRE_INIT_PROB_PARAMS_HPP
#define SWE_PRE_INIT_PROB_PARAMS_HPP

namespace SWE {
void Problem::initialize_problem_parameters(const ProblemInputType& problem_specific_input) {
    SWE::Global::g         = problem_specific_input.g;
    SWE::Global::rho_air   = problem_specific_input.rho_air;
    SWE::Global::rho_water = problem_specific_input.rho_water;

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

    if (problem_specific_input.tidal_potential.type != SWE::TidalPotentialType::None) {
        SWE::SourceTerms::tidal_potential = true;
    }

    if (problem_specific_input.coriolis.type != SWE::CoriolisType::None) {
        SWE::SourceTerms::coriolis = true;
    }

    // specify postprocessin parameters
    if (problem_specific_input.wet_dry.type != SWE::WettingDryingType::None) {
        SWE::PostProcessing::wetting_drying = true;
        SWE::PostProcessing::h_o            = problem_specific_input.wet_dry.h_o;
        SWE::PostProcessing::h_o_threshold  = 1.0e6 * SWE::PostProcessing::h_o * std::numeric_limits<double>::epsilon();
    }

    if (problem_specific_input.slope_limit.type != SWE::SlopeLimitingType::None) {
        SWE::PostProcessing::slope_limiting = true;

        if (problem_specific_input.slope_limit.type == SWE::SlopeLimitingType::CockburnShu) {
            SWE::PostProcessing::M  = problem_specific_input.slope_limit.M;
            SWE::PostProcessing::nu = problem_specific_input.slope_limit.nu;
        }
    }
}
}

#endif