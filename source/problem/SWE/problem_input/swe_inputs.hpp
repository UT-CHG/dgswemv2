#ifndef SWE_INPUTS_HPP
#define SWE_INPUTS_HPP

#include "general_definitions.hpp"
#include "problem/SWE/swe_definitions.hpp"
#include <yaml-cpp/yaml.h>

#include "utilities//file_exists.hpp"

namespace SWE {
// Problem specific preprocessing information containers
struct SphericalProjection {
    SphericalProjectionType type = SphericalProjectionType::None;
    double longitude_o           = 0.0;
    double latitude_o            = 0.0;
    double R                     = 6378200.0;
};

struct InitialConditions {
    InitialConditionsType type = InitialConditionsType::Default;
    double ze_initial          = 0.;
    double qx_initial          = 0.;
    double qy_initial          = 0.;
};

// Problem specific bcis information containers
struct TideInput {
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<double> amplitude;
    std::vector<double> phase;
};

struct FlowInput {
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<double> amplitude;
    std::vector<double> phase;
};

struct LeveeInput {
    double H_barrier;
    double C_subcritical;
    double C_supercritical;
};

// Problem specific forcing terms information containers
struct FunctionSource {
    FunctionSourceType type = FunctionSourceType::None;
};

struct BottomFriction {
    BottomFrictionType type = BottomFrictionType::None;
    double coefficient      = 0.0;
    std::string manning_data_file;
};

struct MeteoForcing {
    MeteoForcingType type = MeteoForcingType::None;
    std::string meteo_data_type;
    std::string raw_meteo_data_file;
    std::string meteo_data_file;
    double frequency;
};

struct TidePotential {
    TidePotentialType type = TidePotentialType::None;
};

struct Coriolis {
    CoriolisType type = CoriolisType::None;
};

// Problem specific postprocessing information containers
struct WettingDrying {
    WettingDryingType type = WettingDryingType::None;
    double h_o             = 0.1;
};

struct SlopeLimiting {
    SlopeLimitingType type = SlopeLimitingType::None;
    std::string slope_limiting_type;
    double M  = 1.0e-8;
    double nu = 1.5;
};

// Problem specific inputs
struct Inputs {
    double g         = 9.81;
    double rho_air   = 1.2250;
    double rho_water = 1000.0;

    SphericalProjection spherical_projection;
    InitialConditions initial_conditions;

    std::map<uint, TideInput> tide_bc_data;
    std::map<uint, FlowInput> flow_bc_data;
    std::map<std::pair<uint, uint>, LeveeInput> levee_is_data;

    FunctionSource function_source;
    BottomFriction bottom_friction;
    MeteoForcing meteo_forcing;
    TidePotential tide_potential;
    Coriolis coriolis;

    WettingDrying wet_dry;
    SlopeLimiting slope_limit;

    Inputs() = default;
    Inputs(YAML::Node& swe_node);

    void read_bcis(const std::string& bcis_file);

    YAML::Node as_yaml_node();
};
}

#endif