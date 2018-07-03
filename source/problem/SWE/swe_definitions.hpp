#ifndef SWE_DEFINITIONS_HPP
#define SWE_DEFINITIONS_HPP

namespace SWE {
namespace Global {
static double g __attribute__((unused))         = 9.81;
static double rho_air __attribute__((unused))   = 1.225;
static double rho_water __attribute__((unused)) = 1000.0;
}

namespace SourceTerms {
static bool function_source __attribute__((unused)) = false;
static bool bottom_friction __attribute__((unused)) = false;
static bool meteo_forcing __attribute__((unused))   = false;
static bool tide_potential __attribute__((unused))  = false;
static bool coriolis __attribute__((unused))        = false;

static double Cf __attribute__((unused)) = 0.0;
}

namespace PostProcessing {
static bool wetting_drying __attribute__((unused)) = false;
static bool slope_limiting __attribute__((unused)) = false;

static double h_o __attribute__((unused))           = 0.1;
static double h_o_threshold __attribute__((unused)) = 1.0e5 * std::numeric_limits<double>::epsilon();

// Cockburn-Shu SL parameters
static double M __attribute__((unused))  = 1.0e-8;
static double nu __attribute__((unused)) = 1.5;
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
