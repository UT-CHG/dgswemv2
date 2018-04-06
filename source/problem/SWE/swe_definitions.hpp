#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

namespace SWE {
namespace Global {
static double g         = 9.81;
static double rho_air   = 1.225;
static double rho_water = 1000.0;

static double Cf           = 0.0;
static double h_o          = 0.01;
static double h_o_treshold = std::numeric_limits<double>::epsilon();
}

namespace SourceTerms {
static bool function_source = false;
static bool bottom_friction = false;
static bool meteo_forcing   = false;
static bool tidal_potential = false;
static bool coriolis        = false;
}

enum BoundaryConditions : uchar { land = 0, tidal = 1, flow = 2, distributed = DISTRIBUTED, internal = INTERNAL };
}

#endif
