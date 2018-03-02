#ifndef SWE_INPUTS_HPP
#define SWE_INPUTS_HPP

#include "../../../general_definitions.hpp"
#include <yaml-cpp/yaml.h>

namespace SWE {
enum class InitialConditionsType {
    Constant,
    Function
};

struct InitialConditions {
    InitialConditionsType type = InitialConditionsType::Constant;
    double ze_initial = 0.;
    double qx_initial = 0.;
    double qy_initial = 0.;
};

enum class FunctionSourceType {
    None,
    Test
};

struct FunctionSource {
    FunctionSourceType type = FunctionSourceType::None;
};

enum class BottomFrictionType {
    None,
    Chezy,
    Manning
};

struct BottomFriction {
    BottomFrictionType type = BottomFrictionType::None;
    double coefficient = 0.0;
    std::string manning_data_file;
};

enum class MeteoForcingType {
    None,
    Test
};

struct MeteoForcing {
    MeteoForcingType type = MeteoForcingType::None;
    std::string meteo_data_file;
    double frequency;
};

enum class TidalPotentialType {
    None,
    Test
};

struct TidalPotential {
    TidalPotentialType type = TidalPotentialType::None;
};

enum class CoriolisType {
    None,
    Test
};

struct Coriolis {
    CoriolisType type = CoriolisType::None;
};

struct Inputs {
    Inputs() = default;
    Inputs(YAML::Node& swe_node);

    YAML::Node as_yaml_node();

    double g = 9.81;
    double h_o = 0.01;

    bool parse_input = false;

    InitialConditions initial_conditions;

    FunctionSource function_source;
    BottomFriction bottom_friction;
    MeteoForcing meteo_forcing;
    TidalPotential tidal_potential;
    Coriolis coriolis;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        ar& g& h_o& initial_conditions& function_source& bottom_friction& meteo_forcing& tidal_potential& coriolis&
            parse_input;
    }
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void serialize(Archive& ar, InitialConditions& ic, unsigned) {
    ar& ic.type& ic.ze_initial& ic.qx_initial& ic.qy_initial;
}

template <typename Archive>
void serialize(Archive& ar, FunctionSource& fs, unsigned) {
    ar& fs.type;
}

template <typename Archive>
void serialize(Archive& ar, BottomFriction& bf, unsigned) {
    ar& bf.type& bf.coefficient& bf.manning_data_file;
}

template <typename Archive>
void serialize(Archive& ar, MeteoForcing& mf, unsigned) {
    ar& mf.type& mf.meteo_data_file& mf.frequency;
}

template <typename Archive>
void serialize(Archive& ar, TidalPotential& tp, unsigned) {
    ar& tp.type;
}

template <typename Archive>
void serialize(Archive& ar, Coriolis& c, unsigned) {
    ar& c.type;
};

#endif
}

#endif