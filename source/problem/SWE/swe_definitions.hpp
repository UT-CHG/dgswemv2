#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

#include "swe_boundary_conditions.hpp"

namespace SWE {
namespace Global {
static constexpr double g = 9.81;
static constexpr double Cf = 0.003;
}

enum BoundaryConditions : uchar {
    land = 0,
    tidal = 1,
    flow = 2,
    distributed = DISTRIBUTED,
    internal = INTERNAL
};
}

#endif
