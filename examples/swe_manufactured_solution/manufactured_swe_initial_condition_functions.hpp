#ifndef SWE_INITIAL_CONDITION_FUNCTIONS_HPP
#define SWE_INITIAL_CONDITION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> ic_q(const double t, const Point<2>& pt) {
    /*constexpr double x1 = -PI;
    constexpr double x2 = PI;
    constexpr double y1 = -PI;
    constexpr double y2 = PI;

    constexpr double Ho = 0.2;
    constexpr double zo = 0.025;
    Utilities::ignore(Ho);

    constexpr double w   = 1;
    constexpr double tau = 0;

    double ic_ze = 2 * zo * cos(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) *
                   cos(w * (t + tau)) / (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double ic_qx = zo * sin(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
                   (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double ic_qy = zo * cos(w * (pt[GlobalCoord::x] - x1)) * sin(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
                   (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));*/

    double ic_ze = exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]) - sin(3 * t));
    double ic_qx = cos(pt[GlobalCoord::x] - 4 * t);
    double ic_qy = sin(pt[GlobalCoord::y] + 4 * t);

    StatVector<double, SWE::n_variables> ic_q{ic_ze, ic_qx, ic_qy};

    return ic_q;
}
}

#endif