#include "swe_inputs.hpp"

namespace SWE {
Inputs::Inputs(YAML::Node& swe_node) {
    if (swe_node["gravity"]) {
        g = swe_node["gravity"].as<double>();
    }

    if (swe_node["h_o"]) {
        h_o = swe_node["h_o"].as<double>();
    }

    const std::string malformatted_sp_warning("  Warning sp conditions are mal-formatted. Using default parameters\n");

    if (YAML::Node sp_node = swe_node["spherical_projection"]) {
        if (sp_node["type"]) {
            std::string sp_str = sp_node["type"].as<std::string>();
            if (sp_str == "None") {
                spherical_projection.type = SphericalProjectionType::None;
            } else if (sp_str == "Enable") {
                if (sp_node["polar_o"] && sp_node["azimuthal_o"] && sp_node["R"]) {
                    spherical_projection.type = SphericalProjectionType::Enable;
                    spherical_projection.azimuthal_o = sp_node["polar_o"].as<double>();
                    spherical_projection.polar_o = sp_node["azimuthal_o"].as<double>();
                    spherical_projection.R = sp_node["R"].as<double>();
                } else {
                    std::cerr << malformatted_sp_warning;
                }
            } else {
                std::cerr << malformatted_sp_warning;
            }
        } else {
            std::cerr << malformatted_sp_warning;
        }
    } else {
        std::cerr << malformatted_sp_warning;
    }

    const std::string malformatted_ic_warning(
        "  Warning initial conditions are mal-formatted. Using default parameters\n");

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
                std::cerr << malformatted_ic_warning;
            }
        }
    } else {
        std::cerr << malformatted_ic_warning;
    }

    const std::string malformatted_fsource_warning(
        "  Warning function source is mal-formatted. Using default parameters\n");

    if (YAML::Node func_source = swe_node["function_source"]) {
        if (func_source["type"]) {
            std::string func_source_str = func_source["type"].as<std::string>();
            if (func_source_str == "None") {
                function_source.type = FunctionSourceType::None;
            } else if (func_source_str == "Enable") {
                function_source.type = FunctionSourceType::Enable;
            } else {
                std::cerr << malformatted_fsource_warning;
            }
        } else {
            std::cerr << malformatted_fsource_warning;
        }
    } else {
        std::cout << malformatted_fsource_warning;
    }

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
            } else if (bf_str == "Manning") {
                if (bf_node["coefficient"] && bf_node["input_file"]) {
                    if (bf_node["coefficient"].as<double>() < 0.) {
                        const std::string err_msg("Error: Chezy friction coefficient must be postive\n");
                        throw std::logic_error(err_msg);
                    }
                    bottom_friction.type = BottomFrictionType::Manning;
                    bottom_friction.coefficient = bf_node["coefficient"].as<double>();
                    bottom_friction.manning_data_file = bf_node["input_file"].as<std::string>();
                } else {
                    std::cerr << malformatted_bf_warning;
                }
            } else {
                std::cerr << malformatted_bf_warning;
            }
        } else {
            std::cerr << malformatted_bf_warning;
        }
    } else {
        std::cout << malformatted_bf_warning;
    }

    const std::string malformatted_meteo_warning(
        "  Warning meteo forcing is mal-formatted. Using default parameters\n");

    if (YAML::Node meteo = swe_node["meteo_forcing"]) {
        if (meteo["type"]) {
            std::string meteo_str = meteo["type"].as<std::string>();
            if (meteo_str == "None") {
                meteo_forcing.type = MeteoForcingType::None;
            } else if (meteo_str == "Test") {
                if (meteo["input_file"] && meteo["frequency"]) {
                    meteo_forcing.type = MeteoForcingType::Test;

                    parse_input = true;
                    meteo_forcing.meteo_data_file = meteo["input_file"].as<std::string>();
                    meteo_forcing.frequency = meteo["frequency"].as<double>();
                } else {
                    std::cerr << malformatted_meteo_warning;
                }
            } else {
                std::cerr << malformatted_meteo_warning;
            }
        } else {
            std::cerr << malformatted_meteo_warning;
        }
    } else {
        std::cout << malformatted_meteo_warning;
    }

    const std::string malformatted_tidal_warning(
        "  Warning tidal potential forcing is mal-formatted. Using default parameters\n");

    if (YAML::Node tidal = swe_node["tidal_potential"]) {
        if (tidal["type"]) {
            std::string tidal_str = tidal["type"].as<std::string>();
            if (tidal_str == "None") {
                tidal_potential.type = TidalPotentialType::None;
            } else if (tidal_str == "Test") {
                tidal_potential.type = TidalPotentialType::Test;
            } else {
                std::cerr << malformatted_tidal_warning;
            }
        } else {
            std::cerr << malformatted_tidal_warning;
        }
    } else {
        std::cout << malformatted_tidal_warning;
    }

    const std::string malformatted_coriolis_warning(
        "  Warning coriolis forcing is mal-formatted. Using default parameters\n");

    if (YAML::Node coriolis_node = swe_node["coriolis"]) {
        if (coriolis_node["type"]) {
            std::string coriolis_str = coriolis_node["type"].as<std::string>();
            if (coriolis_str == "None") {
                coriolis.type = CoriolisType::None;
            } else if (coriolis_str == "Enable") {
                coriolis.type = CoriolisType::Enable;
            } else {
                std::cerr << malformatted_coriolis_warning;
            }
        } else {
            std::cerr << malformatted_coriolis_warning;
        }
    } else {
        std::cout << malformatted_coriolis_warning;
    }
}

YAML::Node Inputs::as_yaml_node() {
    YAML::Node ret;
    ret["name"] = "swe";
    ret["gravity"] = g;
    ret["h_o"] = h_o;

    YAML::Node sp_node;
    switch (spherical_projection.type) {
        case SphericalProjectionType::None:
            sp_node["type"] = "None";
            break;
        case SphericalProjectionType::Enable:
            sp_node["type"] = "Enable";
            sp_node["polar_o"] = spherical_projection.polar_o;
            sp_node["azimuthal_o"] = spherical_projection.azimuthal_o;
            sp_node["R"] = spherical_projection.R;
            break;
    }
    ret["spherical_projection"] = sp_node;

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

    YAML::Node func_src_node;
    switch (function_source.type) {
        case FunctionSourceType::None:
            func_src_node["type"] = "None";
            break;
        case FunctionSourceType::Enable:
            func_src_node["type"] = "Enable";
            break;
    }
    ret["function_source"] = func_src_node;

    YAML::Node bf_node;
    switch (bottom_friction.type) {
        case BottomFrictionType::None:
            bf_node["type"] = "None";
            break;
        case BottomFrictionType::Chezy:
            bf_node["type"] = "Chezy";
            bf_node["coefficient"] = bottom_friction.coefficient;
            break;
        case BottomFrictionType::Manning:
            bf_node["type"] = "Manning";
            bf_node["coefficient"] = bottom_friction.coefficient;
            bf_node["input_file"] = bottom_friction.manning_data_file;
            break;
    }
    ret["bottom_friction"] = bf_node;

    YAML::Node meteo_node;
    switch (meteo_forcing.type) {
        case MeteoForcingType::None:
            meteo_node["type"] = "None";
            break;
        case MeteoForcingType::Test:
            meteo_node["type"] = "Test";
            meteo_node["input_file"] = meteo_forcing.meteo_data_file;
            meteo_node["frequency"] = meteo_forcing.frequency;
            break;
    }
    ret["meteo_forcing"] = meteo_node;

    YAML::Node tide_node;
    switch (tidal_potential.type) {
        case TidalPotentialType::None:
            tide_node["type"] = "None";
            break;
        case TidalPotentialType::Test:
            tide_node["type"] = "Test";
            break;
    }
    ret["tidal_potential"] = tide_node;

    YAML::Node coriolis_node;
    switch (coriolis.type) {
        case CoriolisType::None:
            coriolis_node["type"] = "None";
            break;
        case CoriolisType::Enable:
            coriolis_node["type"] = "Enable";
            break;
    }
    ret["coriolis"] = coriolis_node;

    return ret;
}
}