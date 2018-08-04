#ifndef GN_SOURCE_FUNCTIONS_HPP
#define GN_SOURCE_FUNCTIONS_HPP

namespace GN {
inline double source_ze(const double t, const Point<2>& pt) {
    return 0.0;
}

inline double source_qx(const double t, const Point<2>& pt) {
    return 0.0;
}

inline double source_qy(const double t, const Point<2>& pt) {
    return 0.0;
}

inline StatVector<double, GN::n_variables> source_u(const double t, const Point<2>& pt) {
    double source_ze = 0.0;
    double source_qx = 0.0;
    double source_qy = 0.0;

    StatVector<double, GN::n_variables> source_u{source_ze, source_qx, source_qy};

    return source_u;
}
}

#endif
