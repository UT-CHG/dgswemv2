#ifndef PROBLEM_OMPI_FUNCTIONS_HPP
#define PROBLEM_OMPI_FUNCTIONS_HPP

#ifdef SWE_SUPPORT
#ifdef RKDG_SUPPORT

#include "problem/SWE/discretization_RKDG/kernels_preprocessor/rkdg_swe_pre_ompi.hpp"
#include "problem/SWE/discretization_RKDG/kernels_processor/rkdg_swe_proc_ompi_step.hpp"

#endif
#endif
#ifdef GN_SUPPORT
#ifdef EHDG_SUPPORT

#include "problem/Green-Naghdi/discretization_EHDG/kernels_preprocessor/ehdg_gn_pre_ompi.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/kernels_processor/ehdg_gn_proc_ompi_step.hpp"

#endif
#endif

#endif