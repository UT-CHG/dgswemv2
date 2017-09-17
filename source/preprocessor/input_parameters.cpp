#include "input_parameters.hpp"
#include "ADCIRC_reader/adcirc_format.hpp"

InputParameters::InputParameters(const std::string& input_string) {
    YAML::Node input_file = YAML::LoadFile(input_string);
    // Process Mesh information
    {
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

    polynomial_order = input_file["polynomial_order"].as<uint>();
}

InputParameters::InputParameters(const std::string& input_string, const uint locality_id, const uint submesh_id)
    : InputParameters(input_string) {
    mesh_file_path.insert(mesh_file_path.find_last_of("."),
                          '_' + std::to_string(locality_id) + '_' + std::to_string(submesh_id));
}

void InputParameters::ReadMesh() {
    if (mesh_format == "Adcirc") {
        AdcircFormat adcirc_file(mesh_file_path);
        mesh_data = MeshMetaData(adcirc_file);
    } else if (mesh_format == "Meta") {
        mesh_data = MeshMetaData(mesh_file_path);
    }
}

void InputParameters::WriteTo(const std::string& output_filename) {
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
