#ifndef SWE_INPUTS_HPP
#define SWE_INPUTS_HPP

#include "problem/SWE/swe_definitions.hpp"
#include <yaml-cpp/yaml.h>
#include "utilities/file_exists.hpp"

namespace SWE {
// Problem specific preprocessing information containers
struct SphericalProjection {
    SphericalProjectionType type = SphericalProjectionType::None;
    double longitude_o           = 0.0;
    double latitude_o            = 0.0;
    double R                     = 6378200.0;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type
            & longitude_o
            & latitude_o
            & R;
        // clang-format on
    }
#endif
};

struct InitialConditions {
    InitialConditionsType type = InitialConditionsType::Default;
    double ze_initial          = 0.;
    double qx_initial          = 0.;
    double qy_initial          = 0.;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type
            & ze_initial
            & qx_initial
            & qy_initial;
        // clang-format on
    }
#endif
};

// Problem specific bcis information containers
struct TideInput {
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<double> amplitude;
    std::vector<double> phase;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & frequency
            & forcing_fact
            & equilib_arg
            & amplitude
            & phase;
        // clang-format on
    }
#endif
};

struct FlowInput {
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<double> amplitude;
    std::vector<double> phase;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & frequency
            & forcing_fact
            & equilib_arg
            & amplitude
            & phase;
        // clang-format on
    }
#endif
};

struct LeveeInput {
    double H_barrier;
    double C_subcritical;
    double C_supercritical;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & H_barrier
            & C_subcritical
            & C_supercritical;
        // clang-format on
    }
#endif
};

// Problem specific forcing terms information containers
struct FunctionSource {
    FunctionSourceType type = FunctionSourceType::None;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type;
        // clang-format on
    }
#endif
};

struct BottomFriction {
    BottomFrictionType type = BottomFrictionType::None;
    double coefficient      = 0.0;
    std::string manning_data_file;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type
            & coefficient
            & manning_data_file;
        // clang-format on
    }
#endif
};

struct MeteoForcing {
    MeteoForcingType type = MeteoForcingType::None;
    std::string meteo_data_type;
    std::string raw_meteo_data_file;
    std::string meteo_data_file;
    double frequency;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type
            & meteo_data_file
            & frequency;
        // clang-format on
    }
#endif
};

struct TidePotential {
    TidePotentialType type = TidePotentialType::None;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type;
        // clang-format on
    }
#endif
};

struct Coriolis {
    CoriolisType type = CoriolisType::None;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type;
        // clang-format on
    }
#endif
};

// Problem specific postprocessing information containers
struct WettingDrying {
    WettingDryingType type = WettingDryingType::None;
    double h_o             = 0.1;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type
            & h_o;
        // clang-format on
    }
#endif
};

struct SlopeLimiting {
    SlopeLimitingType type = SlopeLimitingType::None;
    std::string slope_limiting_type;
    double M  = 1.0e-8;
    double nu = 1.5;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & type
            & slope_limiting_type
            & M
            & nu;
        // clang-format on
    }
#endif
};

// Problem specific inputs
struct Inputs {
    std::string name;

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

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & g
            & rho_air
            & rho_water
            & spherical_projection
            & initial_conditions
            & tide_bc_data
            & flow_bc_data
            & levee_is_data
            & function_source
            & bottom_friction
            & meteo_forcing
            & tide_potential
            & coriolis
            & wet_dry
            & slope_limit;
        // clang-format on
    }
#endif
};
}
#endif