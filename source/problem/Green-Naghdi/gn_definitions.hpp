#ifndef GN_DEFINITIONS_HPP
#define GN_DEFINITIONS_HPP

#include "utilities/ignore.hpp"

//#define RKDG_SWE
#define EHDG_SWE
//#define IHDG_SWE

#ifdef RKDG_SWE
namespace SWE_SIM = SWE::RKDG;
#endif

#ifdef EHDG_SWE
namespace SWE_SIM = SWE::EHDG;
#endif

#ifdef IHDG_SWE
namespace SWE_SIM = SWE::IHDG;
#endif

namespace GN {
namespace NDParameters {
static double alpha = 1.0;

const bool ignored_vars = Utilities::ignore(alpha);
}

namespace Global {
static double g         = 9.81;
static double rho_air   = 1.225;
static double rho_water = 1000.0;

const bool ignored_vars = Utilities::ignore(g, rho_air, rho_water);
}

constexpr uint n_dimensions = 2;

//#define D_PROJECTION
#define D_GREENGAUSS
//#define D_LEASTSQUARES

#if defined(D_GREENGAUSS) || defined(D_LEASTSQUARES)
#define D_RECONSTRUCTION
//#define D_RECON_INT
//#define D_RECON_LS
#define D_RECON_AVG
#endif

constexpr uint n_du_terms  = 4;
constexpr uint n_ddu_terms = 8;
enum DU : uint { ux = 0, uy = 1, vx = 2, vy = 3 };
enum DDU : uint { uxx = 0, uxy = 1, uyx = 2, uyy = 3, vxx = 4, vxy = 5, vyx = 6, vyy = 7 };

constexpr uint n_ddbath_terms  = 4;
constexpr uint n_dddbath_terms = 8;
enum DDBath : uint { bxx = 0, bxy = 1, byx = 2, byy = 3 };
enum DDDBath : uint { bxxx = 0, bxxy = 1, bxyx = 2, bxyy = 3, byxx = 4, byxy = 5, byyx = 6, byyy = 7 };

/* These must shadow SWE bc types */
enum BoundaryTypes : uchar { land = 0, tide = 1, flow = 2, function = 3, internal = INTERNAL, levee = INTERNAL + 1 };

namespace EHDG {
constexpr uint n_communications = SWE_SIM::n_communications + 4;
enum CommTypes : uchar {
    dc_global_dof_indx = SWE_SIM::n_communications + 0,
    dbath              = SWE_SIM::n_communications + 1,
    derivatives        = SWE_SIM::n_communications + 2,
    global_mat         = SWE_SIM::n_communications + 3
};
}
}

#endif
