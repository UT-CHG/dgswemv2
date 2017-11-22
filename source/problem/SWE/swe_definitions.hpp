#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

namespace SWE {
namespace Global {
static constexpr double g = 9.81;
static constexpr double Cf = 0.0;
static constexpr double h_o = 0.01;
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
