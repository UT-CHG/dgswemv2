#ifndef GN_DEFINITIONS_HPP
#define GN_DEFINITIONS_HPP

#include "utilities/ignore.hpp"

namespace GN {
namespace NDParameters {
static double alpha     = 1.0;
const bool ignored_vars = Utilities::ignore(alpha);
}

namespace Global {
static double g         = 9.81;
static double rho_air   = 1.225;
static double rho_water = 1000.0;
static double R_earth   = 6378200.0;
const bool ignored_vars = Utilities::ignore(g, rho_air, rho_water, R_earth);
}

namespace SourceTerms {
static bool function_source = false;
static bool bottom_friction = false;
static bool meteo_forcing   = false;
static bool tide_potential  = false;
static bool coriolis        = false;

static double Cf = 0.0;

const bool ignored_vars =
    Utilities::ignore(function_source, bottom_friction, meteo_forcing, tide_potential, coriolis, Cf);
}

constexpr uint n_dimensions  = 2;
constexpr uint n_variables   = 3;
constexpr uint n_auxiliaries = 2;

enum Variables : uint { ze = 0, qx = 1, qy = 2 };

enum Auxiliaries : uint { bath = 0, h = 1 };

enum JacobianVariables : uint {
    ze_ze = 0,
    ze_qx = 1,
    ze_qy = 2,
    qx_ze = 3,
    qx_qx = 4,
    qx_qy = 5,
    qy_ze = 6,
    qy_qx = 7,
    qy_qy = 8
};

constexpr uint n_du_terms  = 4;
constexpr uint n_ddu_terms = 8;

enum FirstDerivatives : uint { ux = 0, uy = 1, vx = 2, vy = 3 };

enum SecondDerivatives : uint { uxx = 0, uxy = 1, uyx = 2, uyy = 3, vxx = 4, vxy = 5, vyx = 6, vyy = 7 };

/* These must shadow SWE bc types */
enum BoundaryTypes : uchar { land = 0, tide = 1, flow = 2, internal = INTERNAL };

enum class SphericalProjectionType { None, Enable };

enum class InitialConditionsType { Default, Constant, Function };

enum class FunctionSourceType { None, Enable };

enum class BottomFrictionType { None, Chezy, Manning };

enum class MeteoForcingType { None, Enable };

enum class TidePotentialType { None, Test };  // not yet implemented

enum class CoriolisType { None, Enable };
}

#endif