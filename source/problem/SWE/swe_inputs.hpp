#ifndef SWE_INPUTS
#define SWE_INPUTS

#include <yaml-cpp/yaml.h>

namespace SWE {

enum class BottomFrictionType {
    None,
    Chezy
};

struct BottomFriction {
    BottomFrictionType type = BottomFrictionType::None;
    double coefficient = 0.0;
};

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

struct Inputs {
    double g = 9.81;

    BottomFriction bottom_friction;

    InitialConditions initial_conditions;

    Inputs(YAML::Node& swe_node);
};
}
#endif