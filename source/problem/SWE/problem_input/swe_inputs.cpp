#include "swe_inputs.hpp"

namespace SWE {
Inputs::Inputs(YAML::Node& swe_node) {
    if ( swe_node["name"] ) {
        this->name = swe_node["name"].as<std::string>();

        //check for acceptable names
        if ( this->name != "rkdg_swe" && this->name != "ehdg_swe" && this->name != "ihdg_swe" ) {
            throw std::runtime_error{"Fatal error: unknown problem name: "+this->name+'\n'};
        }
    }

    if (swe_node["gravity"]) {
        this->g = swe_node["gravity"].as<double>();
    } else {
        std::cerr << "Using default g\n";
    }

    if (swe_node["density_air"]) {
        this->rho_air = swe_node["density_air"].as<double>();
    } else {
        std::cerr << "Using default air density\n";
    }

    if (swe_node["density_water"]) {
        this->rho_water = swe_node["density_water"].as<double>();
    } else {
        std::cerr << "Using default water density\n";
    }

    const std::string malformatted_sp_warning(
        "Warning: spherical projection inputs are mal-formatted. Using default parameters.\n");

    if (YAML::Node sp_node = swe_node["spherical_projection"]) {
        if (sp_node["longitude_o"] && sp_node["latitude_o"] && sp_node["R"]) {
            this->spherical_projection.type = SphericalProjectionType::Enable;

            this->spherical_projection.longitude_o = sp_node["longitude_o"].as<double>();
            this->spherical_projection.latitude_o  = sp_node["latitude_o"].as<double>();
            this->spherical_projection.R           = sp_node["R"].as<double>();
        } else {
            std::cerr << malformatted_sp_warning;
        }
    }

    const std::string malformatted_ic_warning(
        "Warning: initial conditions are mal-formatted. Using default parameters.\n");

    if (YAML::Node ic_node = swe_node["initial_conditions"]) {
        if (ic_node["type"]) {
            std::string ic_str = ic_node["type"].as<std::string>();

            if (ic_str == "Constant") {
                if (ic_node["initial_surface_height"] && ic_node["initial_momentum_x"] &&
                    ic_node["initial_momentum_y"]) {
                    this->initial_conditions.type = InitialConditionsType::Constant;

                    this->initial_conditions.ze_initial = ic_node["initial_surface_height"].as<double>();
                    this->initial_conditions.qx_initial = ic_node["initial_momentum_x"].as<double>();
                    this->initial_conditions.qy_initial = ic_node["initial_momentum_y"].as<double>();
                } else {
                    std::cerr << malformatted_ic_warning;
                }
            } else if (ic_str == "Function") {
                this->initial_conditions.type = InitialConditionsType::Function;
            } else {
                std::cerr << malformatted_ic_warning;
            }
        } else {
            std::cerr << malformatted_ic_warning;
        }
    }

    const std::string malformatted_fsource_warning(
        "Warning: function source is mal-formatted. Using default parameters.\n");

    if (YAML::Node func_source = swe_node["function_source"]) {
        if (!func_source.IsNull()) {
            std::string func_source_str = func_source.as<std::string>();

            if (func_source_str == "Enable") {
                this->function_source.type = FunctionSourceType::Enable;
            } else {
                std::cerr << malformatted_fsource_warning;
            }
        } else {
            std::cerr << malformatted_fsource_warning;
        }
    }

    const std::string malformatted_bf_warning("Warning: bottom friction is mal-formatted. Using default parameters.\n");

    if (YAML::Node bf_node = swe_node["bottom_friction"]) {
        if (bf_node["type"]) {
            std::string bf_str = bf_node["type"].as<std::string>();

            if (bf_str == "Chezy") {
                if (bf_node["coefficient"]) {
                    this->bottom_friction.type = BottomFrictionType::Chezy;

                    this->bottom_friction.coefficient = bf_node["coefficient"].as<double>();

                    if (this->bottom_friction.coefficient < 0.) {
                        throw std::logic_error("Fatal Error: Chezy friction coefficient must be postive!\n");
                    }
                } else {
                    std::cerr << malformatted_bf_warning;
                }
            } else if (bf_str == "Manning") {
                if (bf_node["coefficient"] && bf_node["input_file"]) {
                    this->bottom_friction.type = BottomFrictionType::Manning;

                    this->bottom_friction.coefficient       = bf_node["coefficient"].as<double>();
                    this->bottom_friction.manning_data_file = bf_node["input_file"].as<std::string>();

                    if (this->bottom_friction.coefficient < 0.) {
                        throw std::logic_error("Fatal Error: Chezy friction coefficient must be postive!\n");
                    }
                } else {
                    std::cerr << malformatted_bf_warning;
                }
            } else {
                std::cerr << malformatted_bf_warning;
            }
        } else {
            std::cerr << malformatted_bf_warning;
        }
    }

    const std::string malformatted_meteo_warning(
        "Warning: meteo forcing is mal-formatted. Using default parameters.\n");

    if (YAML::Node meteo = swe_node["meteo_forcing"]) {
        if (meteo["type"] && meteo["raw_input_file"] && meteo["input_file"] && meteo["frequency"]) {
            this->meteo_forcing.type = MeteoForcingType::Enable;

            this->meteo_forcing.meteo_data_type     = meteo["type"].as<std::string>();
            this->meteo_forcing.raw_meteo_data_file = meteo["raw_input_file"].as<std::string>();
            this->meteo_forcing.meteo_data_file     = meteo["input_file"].as<std::string>();
            this->meteo_forcing.frequency           = meteo["frequency"].as<double>();
        } else {
            std::cerr << malformatted_meteo_warning;
        }
    }

    const std::string malformatted_tide_warning(
        "Warning: tide potential forcing is mal-formatted. Using default parameters.\n");

    if (YAML::Node tide = swe_node["tide_potential"]) {
        if (tide["type"]) {
            std::string tide_str = tide["type"].as<std::string>();

            if (tide_str == "Test") {
                this->tide_potential.type = TidePotentialType::Test;
            } else {
                std::cerr << malformatted_tide_warning;
            }
        } else {
            std::cerr << malformatted_tide_warning;
        }
    }

    const std::string malformatted_coriolis_warning(
        "Warning: coriolis forcing is mal-formatted. Using default parameters.\n");

    if (YAML::Node coriolis_node = swe_node["coriolis"]) {
        if (!coriolis_node.IsNull()) {
            std::string coriolis_str = coriolis_node.as<std::string>();

            if (coriolis_str == "Enable") {
                this->coriolis.type = CoriolisType::Enable;
            } else {
                std::cerr << malformatted_coriolis_warning;
            }
        } else {
            std::cerr << malformatted_coriolis_warning;
        }
    }

    const std::string malformatted_wd_warning("Warning: wet-dry is mal-formatted. Using default parameters.\n");

    if (YAML::Node wd_node = swe_node["wetting_drying"]) {
        if (wd_node["h_o"]) {
            this->wet_dry.type = WettingDryingType::Enable;

            this->wet_dry.h_o = wd_node["h_o"].as<double>();
        } else {
            std::cerr << malformatted_wd_warning;
        }
    }

    const std::string malformatted_sl_warning("Warning: slope limiting is mal-formatted. Using default parameters.\n");

    if (YAML::Node sl_node = swe_node["slope_limiting"]) {
        if (sl_node["type"]) {
            std::string sl_string = sl_node["type"].as<std::string>();

            if (sl_string == "Cockburn-Shu") {
                if (sl_node["M"] && sl_node["nu"]) {
                    this->slope_limit.type = SlopeLimitingType::CockburnShu;

                    this->slope_limit.slope_limiting_type = "Cockburn-Shu";
                    this->slope_limit.M                   = sl_node["M"].as<double>();
                    this->slope_limit.nu                  = sl_node["nu"].as<double>();
                } else {
                    std::cerr << malformatted_sl_warning;
                }
            } else {
                std::cerr << malformatted_sl_warning;
            }
        } else {
            std::cerr << malformatted_sl_warning;
        }
    }
}

void Inputs::read_bcis(const std::string& bcis_file) {
    /*if (!Utilities::file_exists(bcis_file)) {
        throw std::logic_error("Fata Error: BC/ISP data file : " + bcis_file + " was not found!\n");
    }*/

    std::ifstream file(bcis_file);

    uint nnodes, btype;

    std::string line;
    std::stringstream stream;

    while (std::getline(file, line)) {
        stream = std::stringstream(line);
        stream >> btype >> nnodes;

        if (btype == SWE::BoundaryTypes::tide) {
            uint ncon, node_ID;
            double frequency, force_fact, eq_argument, amplitude, phase;

            std::getline(file, line);

            stream = std::stringstream(line);
            stream >> ncon;

            for (uint con = 0; con < ncon; ++con) {
                std::getline(file, line);

                stream = std::stringstream(line);
                stream >> frequency >> force_fact >> eq_argument;

                for (uint node = 0; node < nnodes; ++node) {
                    std::getline(file, line);

                    stream = std::stringstream(line);
                    stream >> node_ID >> amplitude >> phase;

                    TideInput& tide = this->tide_bc_data[node_ID];

                    tide.frequency.emplace_back(frequency);
                    tide.forcing_fact.emplace_back(force_fact);
                    tide.equilib_arg.emplace_back(eq_argument);

                    tide.amplitude.emplace_back(amplitude);
                    tide.phase.emplace_back(phase);
                }
            }
        } else if (btype == SWE::BoundaryTypes::flow) {
            uint ncon, node_ID;
            double frequency, force_fact, eq_argument, amplitude, phase;

            std::getline(file, line);

            stream = std::stringstream(line);
            stream >> ncon;

            for (uint con = 0; con < ncon; ++con) {
                std::getline(file, line);

                stream = std::stringstream(line);
                stream >> frequency >> force_fact >> eq_argument;

                for (uint node = 0; node < nnodes; ++node) {
                    std::getline(file, line);

                    stream = std::stringstream(line);
                    stream >> node_ID >> amplitude >> phase;

                    FlowInput& flow = this->flow_bc_data[node_ID];

                    flow.frequency.emplace_back(frequency);
                    flow.forcing_fact.emplace_back(force_fact);
                    flow.equilib_arg.emplace_back(eq_argument);

                    flow.amplitude.emplace_back(amplitude);
                    flow.phase.emplace_back(phase);
                }
            }
        } else if (btype == SWE::BoundaryTypes::levee) {
            uint front_node, back_node;

            for (uint node = 0; node < nnodes; ++node) {
                std::getline(file, line);

                LeveeInput levee;

                stream = std::stringstream(line);
                stream >> front_node >> back_node >> levee.H_barrier >> levee.C_subcritical >> levee.C_supercritical;

                this->levee_is_data[std::pair<uint, uint>{front_node, back_node}] = levee;
            }
        }
    }
}

YAML::Node Inputs::as_yaml_node() {
    YAML::Node ret;

    ret["name"]          = this->name;
    ret["gravity"]       = this->g;
    ret["density_air"]   = this->rho_air;
    ret["density_water"] = this->rho_water;

    YAML::Node sp_node;
    switch (this->spherical_projection.type) {
        case SphericalProjectionType::None:
            break;
        case SphericalProjectionType::Enable:
            sp_node["type"]        = "Enable";
            sp_node["longitude_o"] = this->spherical_projection.longitude_o;
            sp_node["latitude_o"]  = this->spherical_projection.latitude_o;
            sp_node["R"]           = this->spherical_projection.R;

            ret["spherical_projection"] = sp_node;
            break;
    }

    YAML::Node ic_node;
    switch (this->initial_conditions.type) {
        case InitialConditionsType::Default:
            break;
        case InitialConditionsType::Constant:
            ic_node["type"]                   = "Constant";
            ic_node["initial_surface_height"] = this->initial_conditions.ze_initial;
            ic_node["initial_momentum_x"]     = this->initial_conditions.qx_initial;
            ic_node["initial_momentum_y"]     = this->initial_conditions.qy_initial;

            ret["initial_conditions"] = ic_node;
            break;
        case InitialConditionsType::Function:
            ic_node["type"] = "Function";

            ret["initial_conditions"] = ic_node;
            break;
    }

    switch (this->function_source.type) {
        case FunctionSourceType::None:
            break;
        case FunctionSourceType::Enable:
            ret["function_source"] = "Enable";
            break;
    }

    YAML::Node bf_node;
    switch (this->bottom_friction.type) {
        case BottomFrictionType::None:
            break;
        case BottomFrictionType::Chezy:
            bf_node["type"]        = "Chezy";
            bf_node["coefficient"] = this->bottom_friction.coefficient;

            ret["bottom_friction"] = bf_node;
            break;
        case BottomFrictionType::Manning:
            bf_node["type"]        = "Manning";
            bf_node["coefficient"] = this->bottom_friction.coefficient;
            bf_node["input_file"]  = this->bottom_friction.manning_data_file;

            ret["bottom_friction"] = bf_node;
            break;
    }

    YAML::Node meteo_node;
    switch (this->meteo_forcing.type) {
        case MeteoForcingType::None:
            break;
        case MeteoForcingType::Enable:
            meteo_node["type"]           = this->meteo_forcing.meteo_data_type;
            meteo_node["raw_input_file"] = this->meteo_forcing.raw_meteo_data_file;
            meteo_node["input_file"]     = this->meteo_forcing.meteo_data_file;
            meteo_node["frequency"]      = this->meteo_forcing.frequency;

            ret["meteo_forcing"] = meteo_node;
            break;
    }

    YAML::Node tide_node;
    switch (this->tide_potential.type) {
        case TidePotentialType::None:
            break;
        case TidePotentialType::Test:
            tide_node["type"] = "Test";

            ret["tide_potential"] = tide_node;
            break;
    }

    switch (this->coriolis.type) {
        case CoriolisType::None:
            break;
        case CoriolisType::Enable:
            ret["coriolis"] = "Enable";
            break;
    }

    YAML::Node wd_node;
    switch (this->wet_dry.type) {
        case WettingDryingType::None:
            break;
        case WettingDryingType::Enable:
            wd_node["h_o"] = this->wet_dry.h_o;

            ret["wetting_drying"] = wd_node;
            break;
    }

    YAML::Node sl_node;
    switch (this->slope_limit.type) {
        case SlopeLimitingType::None:
            break;
        case SlopeLimitingType::CockburnShu:
            sl_node["type"] = this->slope_limit.slope_limiting_type;
            sl_node["M"]    = this->slope_limit.M;
            sl_node["nu"]   = this->slope_limit.nu;

            ret["slope_limiting"] = sl_node;
            break;
    }

    return ret;
}
}