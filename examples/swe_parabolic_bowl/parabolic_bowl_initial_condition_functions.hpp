#ifndef SWE_INITIAL_CONDITION_FUNCTIONS_HPP
#define SWE_INITIAL_CONDITION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> ic_q(const double t, const Point<2>& pt) {
    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];
    Utilities::ignore(t);

    constexpr double alpha = 1.e-7;
    constexpr double X     = 1;
    constexpr double Y     = -0.41884;
    const double tau       = 2 * PI / (8 * Global::g * alpha);

    double r2                  = x * x + y * y;
    double bathymetry          = -alpha * r2;
    double water_column_height = 1 / (X + Y) + alpha * (Y * Y - X * X) / ((X + Y) * (X + Y)) * r2;

    // Note that the solution is non-zero for
    // r < sqrt( (X+Y)/(alpha*(X*X - Y*Y)))
    // However, the initialization will check for negative water column height
    // and correct it.
    StatVector<double, SWE::n_variables> ic_q{water_column_height - bathymetry, 0, 0};

    return ic_q;
}
}

#endif
