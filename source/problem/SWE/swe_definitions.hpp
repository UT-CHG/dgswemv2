#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

#include "swe_data.hpp"
#include "swe_boundary_conditions.hpp"
#include "../../geometry/mesh_definitions.hpp"

namespace SWE {
namespace Global {
static constexpr double g = 9.81;
static constexpr double Cf = 0;
}

enum BoundaryConditions : uchar {
    land = 0,
    tidal = 1,
    internal = INTERNAL
};
}

#endif
