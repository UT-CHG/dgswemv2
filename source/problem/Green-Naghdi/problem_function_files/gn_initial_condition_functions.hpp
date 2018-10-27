#ifndef GN_INITIAL_CONDITION_FUNCTIONS_HPP
#define GN_INITIAL_CONDITION_FUNCTIONS_HPP

namespace GN {
inline StatVector<double, SWE::n_variables> ic_u(const double t, const Point<2>& pt) {
    constexpr double g = 9.81;

    constexpr double ao = 0.1;
    constexpr double Ho = 0.3;
    constexpr double xo = 0.0;

    constexpr double co = std::sqrt(g * (Ho + ao));
    constexpr double w  = std::sqrt(3.0 * ao) / (2.0 * Ho * std::sqrt(Ho + ao));

    double ic_ze = ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double ic_qx = co * ao * std::pow(1.0 / std::cosh(w * (pt[GlobalCoord::x] - xo - co * t)), 2);

    double ic_qy = 0.0;

    StatVector<double, SWE::n_variables> ic_u{ic_ze, ic_qx, ic_qy};

    return ic_u;
}
}

#endif
