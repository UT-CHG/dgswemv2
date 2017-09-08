#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP

#include "mesh_metadata.hpp"

#include <yaml-cpp/yaml.h>

struct RKInput {
    uint nstages;
    uint order;
};

struct InputParameters {
    InputParameters() = default;
    InputParameters(const std::string& input_string);
    InputParameters(const std::string& input_string, const uint locality_id, const uint submesh_id);

    void ReadMesh();
    void WriteTo(const std::string& output_filename);

    std::string mesh_file_name;
    std::string mesh_format;
    MeshMetaData mesh_data;

    // right now we only support SSPRK timestepping
    RKInput rk;

    // time parameters
    double dt;
    // double T_start;
    double T_end;

    uint polynomial_order;
};

#endif
