#include "swe_inputs.hpp"

namespace SWE {
Inputs::Inputs(YAML::Node& swe_node) {
    if (swe_node["gravity"]) {
        g = swe_node["gravity"].as<double>();
    }
    // std::cout << " Gravity due to acceleration set to " + std::to_string(g) + '\n';

    const std::string malformatted_bf_warning("  Warning bottom friction is mal-formatted. Using default parameters\n");

    if (YAML::Node bf_node = swe_node["bottom_friction"]) {
        if (bf_node["type"]) {
            std::string bf_str = bf_node["type"].as<std::string>();
            if (bf_str == "None") {
                bottom_friction.type = BottomFrictionType::None;
            } else if (bf_str == "Chezy") {
                if (bf_node["coefficient"]) {
                    if (bf_node["coefficient"].as<double>() < 0.) {
                        const std::string err_msg("Error: Chezy friction coefficient must be postive\n");
                        throw std::logic_error(err_msg);
                    }
                    bottom_friction.type = BottomFrictionType::Chezy;
                    bottom_friction.coefficient = bf_node["coefficient"].as<double>();
                } else {
                    std::cerr << malformatted_bf_warning;
                }
            } else {
                std::cerr << "  Unable to determine bottom friction type; using default parameters\n";
            }
        } else {
            std::cerr << malformatted_bf_warning;
        }
    } else {
        std::cout << "  Bottom Friction unset; using default parameters\n";
    }

    const std::string malformatted_ic_warning(
        "  Warning initial conditions are mal-formatted. Using default parameters");

    if (YAML::Node ic_node = swe_node["initial_conditions"]) {
        if (ic_node["type"]) {
            std::string ic_str = ic_node["type"].as<std::string>();
            if (ic_str == "Constant") {
                if (ic_node["initial_surface_height"] && ic_node["initial_momentum_x"] &&
                    ic_node["initial_momentum_y"]) {
                    initial_conditions.type = InitialConditionsType::Constant;
                    initial_conditions.ze_initial = ic_node["initial_surface_height"].as<double>();
                    initial_conditions.qx_initial = ic_node["initial_momentum_x"].as<double>();
                    initial_conditions.qy_initial = ic_node["initial_momentum_y"].as<double>();
                } else {
                    std::cerr << malformatted_ic_warning;
                }
            } else if (ic_str == "Function") {
                initial_conditions.type = InitialConditionsType::Function;
            } else {
                std::cerr << "  Unable to determine initial conditions type; using default parameters\n";
            }
        }
    } else {
        std::cerr << "  Initial conditions unset; using default parameters\n";
    }
}

YAML::Node Inputs::as_yaml_node() {
    YAML::Node ret;
    ret["name"] = "swe";
    ret["gravity"] = g;

    YAML::Node bf_node;
    switch (bottom_friction.type) {
        case BottomFrictionType::None:
            bf_node["type"] = "None";
            break;
        case BottomFrictionType::Chezy:
            bf_node["type"] = "Chezy";
            bf_node["coefficient"] = bottom_friction.coefficient;
            break;
    }
    ret["bottom_friction"] = bf_node;

    YAML::Node ic_node;
    switch (initial_conditions.type) {
        case InitialConditionsType::Constant:
            ic_node["type"] = "Constant";
            ic_node["initial_surface_height"] = initial_conditions.ze_initial;
            ic_node["initial_momentum_x"] = initial_conditions.qx_initial;
            ic_node["initial_momentum_y"] = initial_conditions.qy_initial;
            break;
        case InitialConditionsType::Function:
            ic_node["type"] = "Function";
            break;
    }
    ret["initial_conditions"] = ic_node;

    return ret;
}
}