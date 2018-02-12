#ifndef SWE_TRUE_SOLUTION_FUNCTIONS_HPP
#define SWE_TRUE_SOLUTION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline double true_ze(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;
    Utilities::ignore(Ho);

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    return 2 * zo * cos(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) * cos(w * (t + tau)) /
           (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));
}

inline double true_qx(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;
    Utilities::ignore(Ho);

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    return zo * sin(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
           (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));
}

inline double true_qy(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;
    Utilities::ignore(Ho);

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    return zo * cos(w * (pt[GlobalCoord::x] - x1)) * sin(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
           (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));
}
}

#endif
