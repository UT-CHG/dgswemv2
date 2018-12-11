#ifndef PROBLEM_FUNCTION_FILES_FWD
#define PROBLEM_FUNCTION_FILES_FWD

#include "general_definitions.hpp"
#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
StatVector<double, SWE::n_variables> ic_q(const double t, const Point<2>& pt);

StatVector<double, SWE::n_variables> source_q(const double t, const Point<2>& pt);

StatVector<double, SWE::n_variables> true_q(const double t, const Point<2>& pt);
}

#endif