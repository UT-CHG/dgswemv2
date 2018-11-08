#ifndef PROBLEM_DEFINITIONS_HPP
#define PROBLEM_DEFINITIONS_HPP

namespace SWE {
namespace RKDG {
struct Problem;
}
namespace EHDG {
struct Problem;
}
namespace IHDG {
struct Problem;
}
}

namespace GN {
namespace EHDG {
struct Problem;
}
}

#ifdef SWE_SUPPORT
#ifdef RKDG_SUPPORT
#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"
#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_kernels_preprocessor.hpp"
#endif
#ifdef EHDG_SUPPORT
#include "problem/SWE/discretization_EHDG/ehdg_swe_problem.hpp"
#include "problem/SWE/discretization_EHDG/kernels_preprocessor/ehdg_swe_kernels_preprocessor.hpp"
#endif
#ifdef IHDG_SUPPORT
#include "problem/SWE/discretization_IHDG/ihdg_swe_problem.hpp"
#include "problem/SWE/discretization_IHDG/kernels_preprocessor/ihdg_swe_kernels_preprocessor.hpp"
#endif
#endif

#ifdef GN_SUPPORT
#ifdef EHDG_SUPPORT
#include "problem/Green-Naghdi/discretization_EHDG/ehdg_gn_problem.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/kernels_preprocessor/ehdg_gn_kernels_preprocessor.hpp"
#endif
#endif

#endif