#ifndef SWE_SOURCE_FUNCTIONS_HPP
#define SWE_SOURCE_FUNCTIONS_HPP

#include "general_definitions.hpp"

namespace SWE {
inline double source_ze(const double t, const Point<2>& pt) {
    return 0;
}

inline double source_qx(const double t, const Point<2>& pt) {
    return 0;
}

inline double source_qy(const double t, const Point<2>& pt) {
    return 0;
}

inline Vector<double, SWE::n_variables> source_u(const double t, const Point<2>& pt) {
    double source_ze = 0.0;
    double source_qx = 0.0;
    double source_qy = 0.0;

    Vector<double, SWE::n_variables> source_u{source_ze, source_qx, source_qy};

    return source_u;
}
}

#endif
