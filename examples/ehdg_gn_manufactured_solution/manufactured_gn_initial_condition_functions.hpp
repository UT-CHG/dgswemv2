#ifndef GN_INITIAL_CONDITION_FUNCTIONS_HPP
#define GN_INITIAL_CONDITION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace GN {
inline double ic_ze(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.025;
    constexpr double Ho = 0.5;
    constexpr double xo = -4.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    return ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);
}

inline double ic_qx(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.025;
    constexpr double Ho = 0.5;
    constexpr double xo = -4.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    return co * ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);
}

inline double ic_qy(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.025;
    constexpr double Ho = 0.5;
    constexpr double xo = -4.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    Utilities::ignore(g, ao, Ho, xo, co, w);

    return 0.0;
}

inline StatVector<double, GN::n_variables> ic_u(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.025;
    constexpr double Ho = 0.5;
    constexpr double xo = -4.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    double ic_ze = ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double ic_qx = co * ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double ic_qy = 0.0;

    StatVector<double, GN::n_variables> ic_u{ic_ze, ic_qx, ic_qy};

    return ic_u;
}
}

#endif
