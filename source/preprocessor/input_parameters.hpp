#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP

#include "general_definitions.hpp"
#include <yaml-cpp/yaml.h>

#include "ADCIRC_reader/adcirc_format.hpp"
#include "mesh_metadata.hpp"

#include "problem/SWE/problem_input/swe_inputs.hpp"
#include "problem/Green-Naghdi/problem_input/gn_inputs.hpp"

struct YamlNodeWrapper {
    YAML::Node node;

    YAML::Node as_yaml_node() { return this->node; }
};

struct MeshInput {
    std::string mesh_format;
    std::string mesh_file_name;
    std::string bc_is_file_name;
    std::string db_file_name;

    CoordinateSystem mesh_coordinate_sys;

    MeshMetaData mesh_data;
    DistributedBoundaryMetaData dbmd_data;
};

struct StepperInput {
    struct tm T_start;
    struct tm T_end;

    double run_time;
    double dt;

    uint nstages;
    uint order;

    double ramp_duration;
};

struct WriterInput {
    bool writing_output{false};
    std::string output_path;

    bool writing_log_file{false};
    bool verbose_log_file{false};
    std::string log_file_name;

    bool writing_vtk_output{false};
    double vtk_output_frequency{std::numeric_limits<double>::max()};
    uint vtk_output_freq_step{std::numeric_limits<uint>::max()};

    bool writing_vtu_output{false};
    double vtu_output_frequency{std::numeric_limits<double>::max()};
    uint vtu_output_freq_step{std::numeric_limits<uint>::max()};

    bool writing_modal_output{false};
    double modal_output_frequency{std::numeric_limits<double>::max()};
    uint modal_output_freq_step{std::numeric_limits<uint>::max()};
};

struct LoadBalancerInput {
    bool use_load_balancer{false};
    std::string name;
    double rebalance_frequency;
};

template <typename ProblemInput = YamlNodeWrapper>
struct InputParameters {
    uint polynomial_order;

    MeshInput mesh_input;
    StepperInput stepper_input;
    ProblemInput problem_input;
    WriterInput writer_input;
    LoadBalancerInput load_balancer_input;

    InputParameters() = default;
    InputParameters(const std::string& input_string);
    InputParameters(const std::string& input_string, const uint locality_id, const uint submesh_id);

    void read_mesh();
    void read_bcis();
    void read_dbmd(const uint locality_id, const uint submesh_id);

    void write_to(const std::string& output_filename);

  private:
    ProblemInput problem_specific_ctor_helper(const YAML::Node& input_file);
};

template <typename ProblemInput>
InputParameters<ProblemInput>::InputParameters(const std::string& input_string) {
    YAML::Node input_file = YAML::LoadFile(input_string);

    // Process P input
    if (input_file["polynomial_order"]) {
        this->polynomial_order = input_file["polynomial_order"].as<uint>();
    } else {
        std::string err_msg{"Error: P YAML node not specified\n"};
        throw std::logic_error(err_msg);
    }

    // Process Mesh information
    if (input_file["mesh"]) {
        YAML::Node raw_mesh = input_file["mesh"];

        if (raw_mesh["format"] && raw_mesh["file_name"] && raw_mesh["coordinate_system"]) {
            this->mesh_input.mesh_format = raw_mesh["format"].as<std::string>();

            if (!((this->mesh_input.mesh_format == "Adcirc") || (this->mesh_input.mesh_format == "Meta"))) {
                std::string err_msg = "Error: Unsupported mesh format: " + this->mesh_input.mesh_format + '\n';
                throw std::logic_error(err_msg);
            }

            this->mesh_input.mesh_file_name = raw_mesh["file_name"].as<std::string>();
            this->mesh_input.bc_is_file_name =
                this->mesh_input.mesh_file_name.substr(0, this->mesh_input.mesh_file_name.find_last_of('.')) + ".bcis";
            this->mesh_input.db_file_name =
                this->mesh_input.mesh_file_name.substr(0, this->mesh_input.mesh_file_name.find_last_of('.')) + ".dbmd";
            std::string coord_sys_string = raw_mesh["coordinate_system"].as<std::string>();

            if (coord_sys_string == "cartesian") {
                this->mesh_input.mesh_coordinate_sys = CoordinateSystem::cartesian;
            } else if (coord_sys_string == "polar") {
                this->mesh_input.mesh_coordinate_sys = CoordinateSystem::polar;
            } else if (coord_sys_string == "spherical") {
                this->mesh_input.mesh_coordinate_sys = CoordinateSystem::spherical;
            } else {
                std::string err_msg = "Error: Unsupported coordinate system: " + coord_sys_string + '\n';
                throw std::logic_error(err_msg);
            }
        } else {
            std::string err_msg{"Error: Mesh YAML node is malformatted\n"};
            throw std::logic_error(err_msg);
        }
    } else {
        std::string err_msg{"Error: Mesh YAML node not specified\n"};
        throw std::logic_error(err_msg);
    }

    // Process timestepping information
    if (input_file["timestepping"]) {
        YAML::Node time_stepping = input_file["timestepping"];

        if (time_stepping["start_time"] && time_stepping["end_time"] && time_stepping["dt"] && time_stepping["order"] &&
            time_stepping["nstages"]) {
            std::string start_time = time_stepping["start_time"].as<std::string>();
            std::string end_time   = time_stepping["end_time"].as<std::string>();

            this->stepper_input.T_start = {0};
            this->stepper_input.T_end   = {0};

            strptime(start_time.c_str(), "%d-%m-%Y %H:%M:%S", &this->stepper_input.T_start);
            strptime(end_time.c_str(), "%d-%m-%Y %H:%M:%S", &this->stepper_input.T_end);

            this->stepper_input.run_time =
                difftime(timegm(&this->stepper_input.T_end), timegm(&this->stepper_input.T_start));
            this->stepper_input.dt      = time_stepping["dt"].as<double>();
            this->stepper_input.order   = time_stepping["order"].as<uint>();
            this->stepper_input.nstages = time_stepping["nstages"].as<uint>();

            this->stepper_input.ramp_duration =
                time_stepping["ramp_duration"] ? time_stepping["ramp_duration"].as<double>() : 0;
        } else {
            std::string err_msg{"Error: Timestepping YAML node is malformatted\n"};
            throw std::logic_error(err_msg);
        }
    } else {
        std::string err_msg{"Error: Timestepping YAML node not specified\n"};
        throw std::logic_error(err_msg);
    }

    // Process problem specific input
    if (input_file["problem"]) {
        this->problem_input = problem_specific_ctor_helper(input_file);
    } else {
        std::string err_msg{"Error: Problem YAML node not specified\n"};
        throw std::logic_error(err_msg);
    }

    // Process dynamic load balancing information
    if (input_file["load_balancer"]) {
        YAML::Node lb_node = input_file["load_balancer"];
        if (lb_node["name"] && lb_node["rebalance_frequency"]) {
            this->load_balancer_input.use_load_balancer   = true;
            this->load_balancer_input.name                = lb_node["name"].as<std::string>();
            this->load_balancer_input.rebalance_frequency = lb_node["rebalance_frequency"].as<double>();
        } else {
            std::string err_msg{"Error: Load Balancer YAML node is malformatted\n"};
            throw std::logic_error(err_msg);
        }
    }

    // Process output information (else no output)
    if (input_file["output"]) {
        YAML::Node out_node = input_file["output"];

        if (out_node["path"]) {
            this->writer_input.writing_output = true;
            this->writer_input.output_path    = out_node["path"].as<std::string>();

            if (this->writer_input.output_path.back() != '/') {
                this->writer_input.output_path += "/";
            }
        } else {
            std::string err_msg{"Error: Output YAML node is malformatted\n"};
            throw std::logic_error(err_msg);
        }

        if (out_node["logfile"]) {
            if (out_node["logfile"]["verbose"] && out_node["logfile"]["name"]) {
                this->writer_input.writing_log_file = true;
                this->writer_input.verbose_log_file = out_node["logfile"]["verbose"].as<bool>();
                this->writer_input.log_file_name    = out_node["logfile"]["name"].as<std::string>();
            } else {
                std::string err_msg("Error: Logfile YAML node is malformatted\n");
                throw std::logic_error(err_msg);
            }
        }

        if (out_node["vtk"]) {
            if (out_node["vtk"]["frequency"]) {
                this->writer_input.writing_vtk_output   = true;
                this->writer_input.vtk_output_frequency = out_node["vtk"]["frequency"].as<double>();
                this->writer_input.vtk_output_freq_step =
                    (uint)std::ceil(this->writer_input.vtk_output_frequency / this->stepper_input.dt);
            } else {
                std::string err_msg("Error: VTK YAML node is malformatted\n");
                throw std::logic_error(err_msg);
            }
        }

        if (out_node["vtu"]) {
            if (out_node["vtu"]["frequency"]) {
                this->writer_input.writing_vtu_output   = true;
                this->writer_input.vtu_output_frequency = out_node["vtu"]["frequency"].as<double>();
                this->writer_input.vtu_output_freq_step =
                    (uint)std::ceil(this->writer_input.vtu_output_frequency / this->stepper_input.dt);
            } else {
                std::string err_msg("Error: VTU YAML node is malformatted\n");
                throw std::logic_error(err_msg);
            }
        }

        if (out_node["modal"]) {
            if (out_node["modal"]["frequency"]) {
                this->writer_input.writing_modal_output   = true;
                this->writer_input.modal_output_frequency = out_node["modal"]["frequency"].as<double>();
                this->writer_input.modal_output_freq_step =
                    (uint)std::ceil(this->writer_input.modal_output_frequency / this->stepper_input.dt);
            } else {
                std::string err_msg("Error: Modal YAML node is malformatted\n");
                throw std::logic_error(err_msg);
            }
        }
    }
}

template <typename ProblemInput>
InputParameters<ProblemInput>::InputParameters(const std::string& input_string,
                                               const uint locality_id,
                                               const uint submesh_id)
    : InputParameters(input_string) {
    this->mesh_input.mesh_file_name.insert(this->mesh_input.mesh_file_name.find_last_of("."),
                                           '_' + std::to_string(locality_id) + '_' + std::to_string(submesh_id));

    this->mesh_input.db_file_name.insert(this->mesh_input.db_file_name.find_last_of("."),
                                         '_' + std::to_string(locality_id) + '_' + std::to_string(submesh_id));
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::read_mesh() {
    if (this->mesh_input.mesh_format == "Adcirc") {
        AdcircFormat adcirc_file(this->mesh_input.mesh_file_name);
        this->mesh_input.mesh_data = MeshMetaData(adcirc_file);
    } else if (this->mesh_input.mesh_format == "Meta") {
        this->mesh_input.mesh_data = MeshMetaData(this->mesh_input.mesh_file_name);
    }
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::read_bcis() {
    this->problem_input.read_bcis(this->mesh_input.bc_is_file_name);
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::read_dbmd(const uint locality_id, const uint submesh_id) {
    this->mesh_input.dbmd_data = DistributedBoundaryMetaData(this->mesh_input.db_file_name, locality_id, submesh_id);
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::write_to(const std::string& output_filename) {
    YAML::Emitter output;

    assert(output.good());

    output << YAML::BeginMap;

    output << YAML::Key << "polynomial_order";
    output << YAML::Value << this->polynomial_order;

    // Assemble mesh information
    YAML::Node mesh;
    mesh["format"]    = this->mesh_input.mesh_format;
    mesh["file_name"] = this->mesh_input.mesh_file_name;

    if (this->mesh_input.mesh_coordinate_sys == CoordinateSystem::cartesian) {
        mesh["coordinate_system"] = "cartesian";
    } else if (this->mesh_input.mesh_coordinate_sys == CoordinateSystem::polar) {
        mesh["coordinate_system"] = "polar";
    } else if (this->mesh_input.mesh_coordinate_sys == CoordinateSystem::spherical) {
        mesh["coordinate_system"] = "spherical";
    }

    output << YAML::Key << "mesh";
    output << YAML::Value << mesh;

    // Assemble timestepping information
    std::array<char, 21> start_time_str;
    std::array<char, 21> end_time_str;

    strftime(start_time_str.data(), 21, "%d-%m-%Y %H:%M:%S", &this->stepper_input.T_start);
    strftime(end_time_str.data(), 21, "%d-%m-%Y %H:%M:%S", &this->stepper_input.T_end);

    YAML::Node timestepping;
    timestepping["start_time"]    = std::string(start_time_str.data());
    timestepping["end_time"]      = std::string(end_time_str.data());
    timestepping["dt"]            = this->stepper_input.dt;
    timestepping["order"]         = this->stepper_input.order;
    timestepping["nstages"]       = this->stepper_input.nstages;
    timestepping["ramp_duration"] = this->stepper_input.ramp_duration;

    output << YAML::Key << "timestepping";
    output << YAML::Value << timestepping;

    output << YAML::Key << "problem";
    output << YAML::Value << this->problem_input.as_yaml_node();

    if (this->writer_input.writing_output) {
        YAML::Node writer;

        writer["path"] = this->writer_input.output_path;

        if (this->writer_input.writing_log_file) {
            writer["logfile"]["name"]    = this->writer_input.log_file_name;
            writer["logfile"]["verbose"] = this->writer_input.verbose_log_file;
        }

        if (this->writer_input.writing_vtk_output) {
            writer["vtk"]["frequency"] = this->writer_input.vtk_output_frequency;
        }

        if (this->writer_input.writing_vtu_output) {
            writer["vtu"]["frequency"] = this->writer_input.vtu_output_frequency;
        }

        if (this->writer_input.writing_modal_output) {
            writer["modal"]["frequency"] = this->writer_input.modal_output_frequency;
        }

        output << YAML::Key << "output";
        output << YAML::Value << writer;
    }

    if (this->load_balancer_input.use_load_balancer) {
        YAML::Node load_balancer;
        load_balancer["name"]                = this->load_balancer_input.name;
        load_balancer["rebalance_frequency"] = this->load_balancer_input.rebalance_frequency;

        output << YAML::Key << "load_balancer";
        output << YAML::Value << load_balancer;
    }

    output << YAML::EndMap;

    std::ofstream ofs(output_filename);
    assert(ofs);

    ofs << "###################################################################"
           "############\n"
        << "#\n"
        << "#  DG_HYPER input file\n"
        << "#\n"
        << "###################################################################"
           "############\n\n";

    ofs << output.c_str();
}

template <>
inline YamlNodeWrapper InputParameters<YamlNodeWrapper>::problem_specific_ctor_helper(const YAML::Node& input_file) {
    YamlNodeWrapper temp;

    temp.node = input_file["problem"];

    return temp;
}

template <>
inline SWE::Inputs InputParameters<typename SWE::Inputs>::problem_specific_ctor_helper(const YAML::Node& input_file) {
    if (input_file["problem"]) {
        YAML::Node swe_node = input_file["problem"];

        if (swe_node["name"].as<std::string>() == "rkdg_swe" || swe_node["name"].as<std::string>() == "ehdg_swe" ||
            swe_node["name"].as<std::string>() == "ihdg_swe") {
            return SWE::Inputs(swe_node);
        } else {
            std::string err_msg{"Error: Unknown problem name: " + swe_node["name"].as<std::string>()};
            throw std::logic_error(err_msg);
        }
    } else {
        std::string err_msg{"Error: Shallow water type must specify a SWE yaml-node\n"};
        throw std::logic_error(err_msg);
    }
}

template <>
inline GN::Inputs InputParameters<typename GN::Inputs>::problem_specific_ctor_helper(const YAML::Node& input_file) {
    if (input_file["problem"]) {
        YAML::Node gn_node = input_file["problem"];
        if (gn_node["name"].as<std::string>() == "ehdg_gn" || gn_node["name"].as<std::string>() == "ihdg_gn") {
            return GN::Inputs(gn_node);
        } else {
            std::string err_msg{"Error: Unknown problem name: " + gn_node["name"].as<std::string>()};
            throw std::logic_error(err_msg);
        }
    } else {
        std::string err_msg{"Error: Green-Naghdi type must specify a GN yaml-node\n"};
        throw std::logic_error(err_msg);
    }
}

#endif
