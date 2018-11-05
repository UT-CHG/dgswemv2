#ifndef PROBLEM_DEFINITIONS_HPP
#define PROBLEM_DEFINITIONS_HPP

namespace SWE {
namespace RKDG {
struct Problem;
}
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
#endif

#ifdef GN_SUPPORT

#endif

#endif