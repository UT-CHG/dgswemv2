#ifndef SWE_INITIAL_CONDITION_FUNCTIONS_HPP
#define SWE_INITIAL_CONDITION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> ic_u(const double t, const Point<2>& pt) {
    constexpr double x1 = 400.;
    constexpr double x2 = 1000.;
    constexpr double y1 = 100.;
    constexpr double y2 = 700.;

    constexpr double Ho = 0.2;
    constexpr double zo = 0.025;
    Utilities::ignore(Ho);

    constexpr double w   = 2 * PI / 600.;
    constexpr double tau = 0;

    double ic_ze = 2 * zo * cos(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) *
                   cos(w * (t + tau)) / (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double ic_qx = zo * sin(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
                   (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double ic_qy = zo * cos(w * (pt[GlobalCoord::x] - x1)) * sin(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
                   (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    StatVector<double, SWE::n_variables> ic_u{ic_ze, ic_qx, ic_qy};

    return ic_u;
}
}

#endif
