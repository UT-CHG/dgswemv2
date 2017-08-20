#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP

#include "mesh_metadata.hpp"

struct RKInput {
    uint nstages;
    uint order;
};

struct InputParameters {
    InputParameters() = default;
    InputParameters(const std::string&);
    InputParameters(const std::string&, uint, uint);

    std::string mesh_file_name;
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
