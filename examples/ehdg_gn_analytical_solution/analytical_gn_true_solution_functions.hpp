#ifndef GN_TRUE_SOLUTION_FUNCTIONS_HPP
#define GN_TRUE_SOLUTION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace GN {
inline StatVector<double, GN::n_variables> true_u(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.05;
    constexpr double Ho = 0.5;
    constexpr double xo = -5.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    double true_ze = ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double true_qx = co * ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double true_qy = 0.0;

    StatVector<double, GN::n_variables> true_u{true_ze, true_qx, true_qy};

    return true_u;
}
}

#endif
