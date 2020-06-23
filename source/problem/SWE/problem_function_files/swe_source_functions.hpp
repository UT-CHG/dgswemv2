#ifndef SWE_SOURCE_FUNCTIONS_HPP
#define SWE_SOURCE_FUNCTIONS_HPP

namespace SWE {
inline StatVector<double, SWE::n_variables> source_q(const double t, const Point<2>& pt) {
    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    double source_ze = -3 * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * cos(3 * t) + cos(4 * t + y) + sin(4 * t - x);

    double source_qx =
        (cos(4 * t - x) * cos(4 * t + y)) / (2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y))) - 4 * sin(4 * t - x) +
        (2 * cos(4 * t - x) * sin(4 * t - x)) / (2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y))) +
        3. * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * (2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y))) *
            cos(3 * x) * sin(3 * y) -
        (3 * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * pow(cos(4 * t - x), 2) * cos(3 * x) * sin(3 * y)) /
            pow(2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)), 2) -
        (3 * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * cos(4 * t - x) * cos(3 * y) * sin(3 * x) * sin(4 * t + y)) /
            pow(2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)), 2);

    double source_qy =
        4 * cos(4 * t + y) +
        3. * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * (2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y))) *
            cos(3 * y) * sin(3 * x) +
        (2 * cos(4 * t + y) * sin(4 * t + y)) / (2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y))) +
        (sin(4 * t - x) * sin(4 * t + y)) / (2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y))) -
        (3 * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * cos(4 * t - x) * cos(3 * x) * sin(3 * y) * sin(4 * t + y)) /
            pow(2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)), 2) -
        (3 * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * cos(3 * y) * sin(3 * x) * pow(sin(4 * t + y), 2)) /
            pow(2 + exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)), 2);

    double source_hc = -3 * exp(-sin(3 * t) + sin(3 * x) * sin(3 * y)) * cos(3 * t) + cos(4 * t + y) + sin(4 * t - x);
    
    StatVector<double, SWE::n_variables> source_q{source_ze, source_qx, source_qy, source_hc};

    return source_q;
}
}

#endif
