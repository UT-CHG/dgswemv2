#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
namespace Global {
static double g         = 9.81;
static double rho_air   = 1.225;
static double rho_water = 1000.0;

const bool ignored_vars = Utilities::ignore(g, rho_air, rho_water);
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

namespace PostProcessing {
static bool wetting_drying = false;
static bool slope_limiting = false;

static double h_o           = 0.1;
static double h_o_threshold = 1.0e-6;

// Cockburn-Shu SL parameters
static double M  = 1.0e-8;
static double nu = 1.5;

const bool ignored_vars = Utilities::ignore(wetting_drying, slope_limiting, h_o, h_o_threshold, M, nu);
}

constexpr uint n_dimensions  = 2;
constexpr uint n_variables   = 3;
constexpr uint n_auxiliaries = 3;

enum Variables : uint { ze = 0, qx = 1, qy = 2 };

enum Auxiliaries : uint { bath = 0, h = 1, sp = 2 };

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

enum BoundaryTypes : uchar { land = 0, tide = 1, flow = 2, internal = INTERNAL, levee = INTERNAL + 1 };

namespace RKDG {
constexpr uint n_communications = 3;
enum CommTypes : uchar { baryctr_coord = 0, bound_state = 1, baryctr_state = 2 };
}

namespace EHDG {
constexpr uint n_communications = 1;
enum CommTypes : uchar { bound_state = 0 };
}

namespace IHDG {
constexpr uint n_communications = 1;
enum CommTypes : uchar { global_dof_indx = 0 };
}

enum class SphericalProjectionType { None, Enable };

enum class InitialConditionsType { Default, Constant, Function };

enum class FunctionSourceType { None, Enable };

enum class BottomFrictionType { None, Chezy, Manning };

enum class MeteoForcingType { None, Enable };

enum class TidePotentialType { None, Test };  // not yet implemented

enum class CoriolisType { None, Enable };

enum class WettingDryingType { None, Enable };

enum class SlopeLimitingType { None, CockburnShu };
}

#endif
