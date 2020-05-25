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
};

struct InitialConditions {
    InitialConditionsType type = InitialConditionsType::Default;
    double ze_initial          = 0.;
    double qx_initial          = 0.;
    double qy_initial          = 0.;
    double hc_initial          = 0.;
};

// Problem specific bcis information containers
struct TideNode {
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<double> amplitude;
    std::vector<double> phase;
};

struct TideBoundary {
    std::map<uint, TideNode> tide_nodes;

    bool get_tide_data(const std::vector<uint>& node_ID, std::vector<TideNode>& tide) {
        for (uint node : node_ID) {  // Check if the boundary segment belongs here
            if (this->tide_nodes.find(node) == this->tide_nodes.end()) {
                return false;
            }
        }

        for (uint node : node_ID) {
            tide.push_back(this->tide_nodes[node]);
        }

        return true;
    }
};

struct FlowNode {
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<double> amplitude;
    std::vector<double> phase;
};

struct FlowBoundary {
    std::map<uint, FlowNode> flow_nodes;

    bool get_flow_data(const std::vector<uint>& node_ID, std::vector<FlowNode>& flow) {
        for (uint node : node_ID) {  // Check if the boundary segment belongs here
            if (this->flow_nodes.find(node) == this->flow_nodes.end()) {
                return false;
            }
        }

        for (uint node : node_ID) {
            flow.push_back(this->flow_nodes[node]);
        }

        return true;
    }
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
    double manning_n        = 0.0;
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

struct SedTransport {
    bool bed_update       = false;
    bool bath_slope_limit = false;
    bool bed_load         = false;
    bool suspended_load   = false;

    double A = 0.0;

    double d              = 0.0;
    double nu             = 0.0;
    double phi            = 0.0;
    double theta_c        = 0.0;
    double rho_sediment   = 0.0;
    double saturation_bed = 0.0;
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
    std::string name;

    double g         = 9.81;
    double rho_air   = 1.2250;
    double rho_water = 1000.0;

    SphericalProjection spherical_projection;
    InitialConditions initial_conditions;

    std::vector<TideBoundary> tide_bc_data;
    std::vector<FlowBoundary> flow_bc_data;
    std::map<std::pair<uint, uint>, LeveeInput> levee_is_data;

    FunctionSource function_source;
    BottomFriction bottom_friction;
    MeteoForcing meteo_forcing;
    TidePotential tide_potential;
    Coriolis coriolis;
    SedTransport sediment_transport;

    WettingDrying wet_dry;
    SlopeLimiting slope_limit;

    Inputs() = default;
    Inputs(YAML::Node& swe_node);

    void read_bcis(const std::string& bcis_file);

    YAML::Node as_yaml_node();
};
}
#endif