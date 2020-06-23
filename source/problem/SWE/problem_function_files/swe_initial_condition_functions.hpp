#ifndef SWE_INITIAL_CONDITION_FUNCTIONS_HPP
#define SWE_INITIAL_CONDITION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> ic_q(const double t, const Point<2>& pt) {
    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    double ic_ze = exp(sin(3 * x) * sin(3 * y) - sin(3 * t));
    double ic_qx = cos(x - 4 * t);
    double ic_qy = sin(y + 4 * t);
    double ic_hc = exp(sin(3 * x) * sin(3 * y) - sin(3 * t)) + 2.0;

    StatVector<double, SWE::n_variables> ic_q{ic_ze, ic_qx, ic_qy, ic_hc};

    return ic_q;
}
}

#endif
