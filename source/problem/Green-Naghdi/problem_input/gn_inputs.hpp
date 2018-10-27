#ifndef GN_INPUTS_HPP
#define GN_INPUTS_HPP

#include "problem/Green-Naghdi/gn_definitions.hpp"
#include <yaml-cpp/yaml.h>
#include "utilities/file_exists.hpp"

namespace GN {
// Problem specific inputs
struct Inputs : SWE::Inputs {
    Inputs() = default;
    Inputs(YAML::Node& gn_node) : SWE::Inputs(gn_node) {}
};
}

#endif