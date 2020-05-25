#ifndef SWE_PRE_INIT_PROB_PARAMS_HPP
#define SWE_PRE_INIT_PROB_PARAMS_HPP

namespace SWE {
void initialize_problem_parameters(const SWE::Inputs& problem_specific_input) {
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

    if (problem_specific_input.tide_potential.type != SWE::TidePotentialType::None) {
        SWE::SourceTerms::tide_potential = true;
    }

    if (problem_specific_input.coriolis.type != SWE::CoriolisType::None) {
        SWE::SourceTerms::coriolis = true;
    }

    if (problem_specific_input.sediment_transport.bed_update) {
        SWE::SedimentTransport::bed_update = true;
        if (problem_specific_input.sediment_transport.bed_load) {
            SWE::SedimentTransport::bed_load = true;

            SWE::SedimentTransport::A = problem_specific_input.sediment_transport.A;
        }
        if (problem_specific_input.sediment_transport.suspended_load) {
            SWE::SedimentTransport::suspended_load = true;

            SWE::Global::rho_sediment = problem_specific_input.sediment_transport.rho_sediment;
            SWE::Global::sat_sediment = problem_specific_input.sediment_transport.saturation_bed;
            SWE::Global::rho_bed      = (1 - SWE::Global::sat_sediment) * SWE::Global::rho_sediment +
                                   SWE::Global::sat_sediment * SWE::Global::rho_water;

            SWE::SedimentTransport::d       = problem_specific_input.sediment_transport.d;
            SWE::SedimentTransport::v       = problem_specific_input.sediment_transport.nu;
            SWE::SedimentTransport::phi     = problem_specific_input.sediment_transport.phi;
            SWE::SedimentTransport::theta_c = problem_specific_input.sediment_transport.theta_c;
        }
        if (problem_specific_input.sediment_transport.bath_slope_limit) {
            SWE::SedimentTransport::bed_slope_limiting = true;
        }
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