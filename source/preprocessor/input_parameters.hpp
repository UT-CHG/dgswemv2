#ifndef INPUT_HPP
#define INPUT_HPP

#include "mesh_metadata.hpp"

struct RKInput {
    uint nstages;
    uint order;
};

struct InputParameters {
    InputParameters(const char* input_string);

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
