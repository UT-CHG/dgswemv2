#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

namespace SWE {
namespace Global {
static double g         = 9.81;
static double rho_air   = 1.225;
static double rho_water = 1000.0;
}

namespace SourceTerms {
static bool function_source = false;
static bool bottom_friction = false;
static bool meteo_forcing   = false;
static bool tide_potential  = false;
static bool coriolis        = false;

static double Cf = 0.0;
}

namespace PostProcessing {
static bool wetting_drying = false;
static bool slope_limiting = false;

static double h_o           = 0.1;
static double h_o_threshold = 1.0e5 * std::numeric_limits<double>::epsilon();

// Cockburn-Shu SL parameters
static double M  = 1.0e-8;
static double nu = 1.5;
}

enum Variables : uint { ze = 0, qx = 1, qy = 2 };

enum BoundaryTypes : uchar {
    land     = 0,
    tide     = 1,
    flow     = 2,
    internal = INTERNAL,
    levee    = INTERNAL + 1,
    periodic = INTERNAL + 2
};

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
