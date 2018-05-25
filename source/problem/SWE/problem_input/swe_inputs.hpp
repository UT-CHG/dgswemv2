#ifndef SWE_INPUTS_HPP
#define SWE_INPUTS_HPP

#include "../../../general_definitions.hpp"
#include "../swe_definitions.hpp"
#include <yaml-cpp/yaml.h>

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

struct TidalPotential {
    TidalPotentialType type = TidalPotentialType::None;
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

// Problem specific boundary condition information containers
struct Levee {
    std::vector<uint> front_nodes;
    std::vector<uint> back_nodes;
    std::vector<double> barrier_height;
    std::vector<double> C_subcritical;
    std::vector<double> C_supercritical;
};

// Problem specific inputs
struct Inputs {
    double g         = 9.81;
    double rho_air   = 1.2250;
    double rho_water = 1000.0;

    SphericalProjection spherical_projection;
    InitialConditions initial_conditions;

    FunctionSource function_source;
    BottomFriction bottom_friction;
    MeteoForcing meteo_forcing;
    TidalPotential tidal_potential;
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