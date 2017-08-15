#include "input_parameters.hpp"
#include "ADCIRC_reader/adcirc_format.hpp"

#include <yaml-cpp/yaml.h>

InputParameters::InputParameters(const char* input_string) {
    YAML::Node input_file = YAML::LoadFile(input_string);

    // Process Mesh information
    {
        YAML::Node raw_mesh = input_file["mesh"];
        std::string format = raw_mesh["format"].as<std::string>();
        if (format == "Adcirc") {
            mesh_file_name = raw_mesh["file_name"].as<std::string>();
            AdcircFormat adcirc_file(mesh_file_name);
            mesh_data = MeshMetaData(adcirc_file);
        } else {
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
    if (polynomial_order > 10) {
        std::string err_msg =
            "Error: Invalid polynomial order: " + std::to_string(polynomial_order) + " can be at most 10\n";

        throw std::logic_error(err_msg);
    }
}

InputParameters::InputParameters(const char* input_string, uint locality, uint thread) {
    YAML::Node input_file = YAML::LoadFile(input_string);

    // Process Mesh information
    {
        YAML::Node raw_mesh = input_file["mesh"];
        std::string format = raw_mesh["format"].as<std::string>();
        if (format == "Adcirc") {
            mesh_file_name = raw_mesh["file_name"].as<std::string>() + std::to_string(locality) + std::to_string(thread) + ".14";
            
            hpx::cout << mesh_file_name;
            
            AdcircFormat adcirc_file(mesh_file_name);
            mesh_data = MeshMetaData(adcirc_file);
        } else {
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
    if (polynomial_order > 10) {
        std::string err_msg =
            "Error: Invalid polynomial order: " + std::to_string(polynomial_order) + " can be at most 10\n";

        throw std::logic_error(err_msg);
    }
}