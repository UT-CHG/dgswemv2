#ifndef SWE_INITIAL_CONDITION_FUNCTIONS_HPP
#define SWE_INITIAL_CONDITION_FUNCTIONS_HPP

#include "general_definitions.hpp"

namespace SWE {
inline double ic_ze(const double t, const Point<2>& pt) {
    return 0;
}

inline double ic_qx(const double t, const Point<2>& pt) {
    return 0;
}

inline double ic_qy(const double t, const Point<2>& pt) {
    return 0;
}

inline Vector<double, SWE::n_variables> ic_u(const double t, const Point<2>& pt) {
    double ic_ze = 0.0;
    double ic_qx = 0.0;
    double ic_qy = 0.0;

    Vector<double, SWE::n_variables> ic_u{ic_ze, ic_qx, ic_qy};

    return ic_u;
}
}

#endif
