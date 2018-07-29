#ifndef GN_TRUE_SOLUTION_FUNCTIONS_HPP
#define GN_TRUE_SOLUTION_FUNCTIONS_HPP

#include "general_definitions.hpp"

namespace GN {
inline double true_ze(const double t, const Point<2>& pt) {
    return 0;
}

inline double true_qx(const double t, const Point<2>& pt) {
    return 0;
}

inline double true_qy(const double t, const Point<2>& pt) {
    return 0;
}

inline StatVector<double, GN::n_variables> true_u(const double t, const Point<2>& pt) {
    double true_ze = 0.0;
    double true_qx = 0.0;
    double true_qy = 0.0;

    StatVector<double, GN::n_variables> true_u{true_ze, true_qx, true_qy};

    return true_u;
}
}

#endif
