#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP

#include "mesh_metadata.hpp"
#include "ADCIRC_reader/adcirc_format.hpp"

#include <yaml-cpp/yaml.h>

#include "problem/SWE/swe_inputs.hpp"

struct YamlNodeWrapper {
    YAML::Node node;

    YAML::Node as_yaml_node() { return node; }
};

struct RKInput {
    uint nstages;
    uint order;
};

struct WriterInput {
    bool write_output = false;
    std::string output_path;
    bool writing_vtk = false;
    uint vtk_output_frequency;
    bool writing_modal_output = false;
    uint modal_output_frequency;
};

template <typename ProblemInput = YamlNodeWrapper>
struct InputParameters {
    InputParameters() = default;
    InputParameters(const std::string& input_string);
    InputParameters(const std::string& input_string, const uint locality_id, const uint submesh_id);

    void ReadMesh();
    void WriteTo(const std::string& output_filename);

    std::string mesh_file_name;
    std::string mesh_file_path;
    std::string mesh_format;
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

  private:
    ProblemInput problem_specific_ctor_helper(const YAML::Node& input_file);
};

template <typename ProblemInput>
InputParameters<ProblemInput>::InputParameters(const std::string& input_string) {
    YAML::Node input_file = YAML::LoadFile(input_string);
    // Process Mesh information
    if (input_file["mesh"]) {
        YAML::Node raw_mesh = input_file["mesh"];
        mesh_format = raw_mesh["format"].as<std::string>();
        mesh_file_name = raw_mesh["file_name"].as<std::string>();

        std::string path_to_input = input_string;
        path_to_input = path_to_input.substr(0, path_to_input.find_last_of("/\\") + 1);
        mesh_file_path = path_to_input + mesh_file_name;

        if (!((mesh_format == "Adcirc") || (mesh_format == "Meta"))) {
            std::string err_msg = "Error: Unsupported mesh format: " + raw_mesh["format"].as<std::string>() + '\n';
            throw std::logic_error(err_msg);
        }
    } else {
        std::string err_msg{"Error: Mesh YAML node not specified\n"};
        throw std::logic_error(err_msg);
    }

    // Process timestepping information
    {
        YAML::Node time_stepping = input_file["timestepping"];
        dt = time_stepping["dt"].as<double>();
        // T_start    = time_stepping["start_time"].as<double>();
        T_end = time_stepping["end_time"].as<double>();
        rk.nstages = time_stepping["order"].as<uint>();
        rk.order = time_stepping["nstages"].as<uint>();
    }

    // Process output information
    if (input_file["output"]) {
        writer_input.write_output = true;
        YAML::Node out_node = input_file["output"];

        writer_input.output_path = out_node["path"].as<std::string>();
        if (writer_input.output_path.back() != '/') {
            writer_input.output_path += "/";
        }

        if (out_node["vtk"]) {
            writer_input.writing_vtk = true;
            writer_input.vtk_output_frequency = out_node["vtk"]["frequency"].as<uint>();
        }

        if (out_node["modal"]) {
            writer_input.writing_modal_output = true;
            writer_input.modal_output_frequency = out_node["modal"]["frequency"].as<uint>();
        }
    }

    if (!input_file["problem"]) {
        std::string err_msg("Error: Problem node not found\n");
        throw std::logic_error(err_msg);
    }
    problem_input = problem_specific_ctor_helper(input_file);

    polynomial_order = input_file["polynomial_order"].as<uint>();
}

template <typename ProblemInput>
InputParameters<ProblemInput>::InputParameters(const std::string& input_string,
                                               const uint locality_id,
                                               const uint submesh_id)
    : InputParameters(input_string) {
    mesh_file_path.insert(mesh_file_path.find_last_of("."),
                          '_' + std::to_string(locality_id) + '_' + std::to_string(submesh_id));
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::ReadMesh() {
    if (mesh_format == "Adcirc") {
        AdcircFormat adcirc_file(mesh_file_path);
        mesh_data = MeshMetaData(adcirc_file);
    } else if (mesh_format == "Meta") {
        mesh_data = MeshMetaData(mesh_file_path);
    }
}

template <typename ProblemInput>
void InputParameters<ProblemInput>::WriteTo(const std::string& output_filename) {
    YAML::Emitter output;
    assert(output.good());
    output << YAML::BeginMap;
    // Assemble mesh information
    {
        YAML::Node mesh;
        mesh["format"] = mesh_format;
        mesh["file_name"] = mesh_file_name;

        output << YAML::Key << "mesh";
        output << YAML::Value << mesh;
    }
    // Assemble timestepping information
    {
        YAML::Node timestepping;
        timestepping["dt"] = dt;
        timestepping["end_time"] = T_end;
        timestepping["order"] = rk.order;
        timestepping["nstages"] = rk.nstages;

        output << YAML::Key << "timestepping" << YAML::Value << timestepping;
    }

    output << YAML::Key << "polynomial_order" << YAML::Value << polynomial_order;

    output << YAML::Key << "problem" << YAML::Value << problem_input.as_yaml_node();

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

// template <typename ProblemInput>
// ProblemInput InputParameter<ProblemInput>::problem_specific_ctor_helper(const YAML::Node& input_file) {
//    static_assert(false, "Error: It seems that for the specific ProblemInput type there is no specialization
// written");
//}

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
