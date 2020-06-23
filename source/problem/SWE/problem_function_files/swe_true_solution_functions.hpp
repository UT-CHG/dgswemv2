#ifndef SWE_TRUE_SOLUTION_FUNCTIONS_HPP
#define SWE_TRUE_SOLUTION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> true_q(const double t, const Point<2>& pt) {
    double true_ze = 0.0;
    double true_qx = 0.0;
    double true_qy = 0.0;
    double true_hc = 0.0;

    StatVector<double, SWE::n_variables> true_q{true_ze, true_qx, true_qy, true_hc};

    return true_q;
}
}

#endif
