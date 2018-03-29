#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP

#include "mesh_metadata.hpp"
#include "ADCIRC_reader/adcirc_format.hpp"

#include <yaml-cpp/yaml.h>

#include "../problem/SWE/problem_input/swe_inputs.hpp"

struct YamlNodeWrapper {
    YAML::Node node;

    YAML::Node as_yaml_node() { return this->node; }
};

struct RKInput {
    uint nstages;
    uint order;
};

struct WriterInput {
    bool        writing_output{false};
    std::string output_path;

    bool        writing_log_file{false};
    bool        verbose_log_file{false};
    std::string log_file_name;

    bool writing_vtk_output{false};
    uint vtk_output_frequency{std::numeric_limits<uint>::max()};

    bool writing_modal_output{false};
    uint modal_output_frequency{std::numeric_limits<uint>::max()};

    YAML::Node as_yaml_node();
};

inline YAML::Node WriterInput::as_yaml_node() {
    YAML::Node ret;

    if (this->writing_output) {
        ret["path"] = this->output_path;
        if (this->writing_log_file) {
            ret["logfile"]["name"]    = this->log_file_name;
            ret["logfile"]["verbose"] = this->verbose_log_file;
        }
        if (this->writing_vtk_output) {
            ret["vtk"]["frequency"] = this->vtk_output_frequency;
        }
        if (this->writing_modal_output) {
            ret["modal"]["frequency"] = this->modal_output_frequency;
        }
    }

    return ret;
}

// Lower case letter for member functions names
template <typename ProblemInput = YamlNodeWrapper>
struct InputParameters {
    std::string  mesh_file_name;
    std::string  mesh_format;
    MeshMetaData mesh_data;

    // right now we only support SSPRK timestepping
    RKInput rk;

    WriterInput writer_input;

    // time parameters
    double dt;
    // double T_start;
    double T_end;

    uint polynomial_order;

    ProblemInput problem_input;

    InputParameters() = default;
    InputParameters(const std::string& input_string);
    InputParameters(const std::string& input_string, const uint locality_id, const uint submesh_id);

    void read_mesh();
    void write_to(const std::string& output_filename);

  private:
    ProblemInput problem_specific_ctor_helper(const YAML::Node& input_file);
};

template <typename ProblemInput>
InputParameters<ProblemInput>::InputParameters(const std::string& input_string) {
    YAML::Node input_file = YAML::LoadFile(input_string);

    // Process Mesh information
    if (input_file["mesh"]) {
        YAML::Node raw_mesh  = input_file["mesh"];
        this->mesh_format    = raw_mesh["format"].as<std::string>();
        this->mesh_file_name = raw_mesh["file_name"].as<std::string>();

        if (!((this->mesh_format == "Adcirc") || (this->mesh_format == "Meta"))) {
            std::string err_msg = "Error: Unsupported mesh format: " + raw_mesh["format"].as<std::string>() + '\n';
            throw std::logic_error(err_msg);
        }
    } else {
        std::string err_msg{"Error: Mesh YAML node not specified\n"};
        throw std::logic_error(err_msg);
    }

    // Process timestepping information
    YAML::Node time_stepping = input_file["timestepping"];
    this->dt                 = time_stepping["dt"].as<double>();
    // T_start    = time_stepping["start_time"].as<double>();
    this->T_end      = time_stepping["end_time"].as<double>();
    this->rk.nstages = time_stepping["order"].as<uint>();
    this->rk.order   = time_stepping["nstages"].as<uint>();

    // Process output information
    if (input_file["output"]) {
        YAML::Node out_node = input_file["output"];

        this->writer_input.writing_output = true;
        this->writer_input.output_path    = out_node["path"].as<std::string>();

        if (this->writer_input.output_path.back() != '/') {
            this->writer_input.output_path += "/";
        }

        if (out_node["logfile"]) {
            this->writer_input.writing_log_file = true;
            this->writer_input.verbose_log_file = out_node["logfile"]["verbose"].as<bool>();
            this->writer_input.log_file_name    = out_node["logfile"]["name"].as<std::string>();
        }

        if (out_node["vtk"]) {
            this->writer_input.writing_vtk_output   = true;
            this->writer_input.vtk_output_frequency = out_node["vtk"]["frequency"].as<uint>();
        }

        if (out_node["modal"]) {
            this->writer_input.writing_modal_output   = true;
            this->writer_input.modal_output_frequency = out_node["modal"]["frequency"].as<uint>();
        }
    }

    if (!input_file["problem"]) {
        std::string err_msg("Error: Problem node not found\n");
        throw std::logic_error(err_msg);
    }

    this->problem_input = problem_specific_ctor_helper(input_file);

    this->polynomial_order = input_file["polynomial_order"].as<uint>();
}

template <typename ProblemInput>
InputParameters<ProblemInput>::InputParameters(const std::string& input_string,
                                               const uint         locality_id,
                                               const uint         submesh_id)
    : InputParameters(input_string) {
    this->mesh_file_name.insert(this->mesh_file_name.find_last_of("."),
                                '_' + std::to_string(locality_id) + '_' + std::to_string(submesh_id));
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::read_mesh() {
    if (this->mesh_format == "Adcirc") {
        AdcircFormat adcirc_file(this->mesh_file_name);
        this->mesh_data = MeshMetaData(adcirc_file);
    } else if (this->mesh_format == "Meta") {
        this->mesh_data = MeshMetaData(this->mesh_file_name);
    }
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::write_to(const std::string& output_filename) {
    YAML::Emitter output;

    assert(output.good());

    output << YAML::BeginMap;

    // Assemble mesh information
    YAML::Node mesh;
    mesh["format"]    = this->mesh_format;
    mesh["file_name"] = this->mesh_file_name;

    output << YAML::Key << "mesh";
    output << YAML::Value << mesh;

    // Assemble timestepping information
    YAML::Node timestepping;
    timestepping["dt"]       = this->dt;
    timestepping["end_time"] = this->T_end;
    timestepping["order"]    = this->rk.order;
    timestepping["nstages"]  = this->rk.nstages;

    output << YAML::Key << "timestepping" << YAML::Value << timestepping;

    output << YAML::Key << "polynomial_order" << YAML::Value << this->polynomial_order;

    output << YAML::Key << "problem" << YAML::Value << this->problem_input.as_yaml_node();

    if (this->writer_input.writing_output) {
        output << YAML::Key << "output" << YAML::Value << this->writer_input.as_yaml_node();
    }

    output << YAML::EndMap;

    std::ofstream ofs(output_filename);
    assert(ofs);

    ofs << "###################################################################"
           "############\n"
        << "#\n"
        << "#  DGSWEMv2 input file\n"
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
    YAML::Node swe_node = input_file["problem"];
    if (swe_node["name"].as<std::string>() == "swe") {
        return SWE::Inputs(swe_node);
    } else {
        std::string err_msg("Error: Shallow water type must specify a SWE yaml-node\n");
        throw std::logic_error(err_msg);
    }
}

#endif