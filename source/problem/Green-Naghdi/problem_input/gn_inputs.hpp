#ifndef GN_INPUTS_HPP
#define GN_INPUTS_HPP

#include "general_definitions.hpp"
#include "problem/Green-Naghdi/gn_definitions.hpp"
#include <yaml-cpp/yaml.h>

#include "utilities/file_exists.hpp"

namespace GN {
// Problem specific inputs
struct Inputs : SWE::Inputs {
    Inputs() = default;
    Inputs(YAML::Node& gn_node) : SWE::Inputs(gn_node) {}

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & g
            & rho_air
            & rho_water
            & initial_conditions
            & tide_bc_data
            & flow_bc_data
            & function_source
            & bottom_friction
            & meteo_forcing
            & tide_potential
            & coriolis;
        // clang-format on
    }
#endif
};
}

#endif