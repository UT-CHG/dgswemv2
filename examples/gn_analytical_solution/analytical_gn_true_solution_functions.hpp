#ifndef GN_TRUE_SOLUTION_FUNCTIONS_HPP
#define GN_TRUE_SOLUTION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> true_q(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.05;
    constexpr double Ho = 0.5;
    constexpr double xo = -5.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    double true_ze = ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double true_qx = co * ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double true_qy = 0.0;

    double true_hc = 0.0;

    StatVector<double, SWE::n_variables> true_q{true_ze, true_qx, true_qy, true_hc};

    return true_q;
}
}

#endif
