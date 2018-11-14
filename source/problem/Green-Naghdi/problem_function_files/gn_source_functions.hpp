#ifndef GN_SOURCE_FUNCTIONS_HPP
#define GN_SOURCE_FUNCTIONS_HPP

namespace GN {
inline StatVector<double, SWE::n_variables> source_q(const double t, const Point<2>& pt) {
    double source_ze = 0.0;
    double source_qx = 0.0;
    double source_qy = 0.0;

    StatVector<double, SWE::n_variables> source_q{source_ze, source_qx, source_qy};

    return source_q;
}
}

#endif
