#ifndef PROBLEM_DEFINITIONS_HPP
#define PROBLEM_DEFINITIONS_HPP

#ifdef SWE_SUPPORT
#ifdef RKDG_SUPPORT

#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"
#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_kernels_preprocessor.hpp"
#endif
#endif

#ifdef GN_SUPPORT

#endif

#endif