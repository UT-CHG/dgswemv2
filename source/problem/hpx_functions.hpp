#ifndef PROBLEM_HPX_FUNCTIONS_HPP
#define PROBLEM_HPX_FUNCTIONS_HPP

#ifdef SWE_SUPPORT
#ifdef RKDG_SUPPORT

#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_pre_hpx.hpp"
#include "problem/SWE/discretization_RKDG/kernels_processor/rkdg_swe_proc_hpx_stage.hpp"

#endif
#ifdef EHDG_SUPPORT

#include "problem/SWE/discretization_EHDG/kernels_preprocessor/ehdg_swe_pre_hpx.hpp"
#include "problem/SWE/discretization_EHDG/kernels_processor/ehdg_swe_proc_hpx_stage.hpp"

#endif
#endif

#endif