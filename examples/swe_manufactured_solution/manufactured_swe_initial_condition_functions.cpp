#include "utilities/ignore.hpp"
#include "problem/SWE/problem_function_files/fwd.hpp"

namespace SWE {
StatVector<double, SWE::n_variables> ic_q(const double t, const Point<2>& pt) {
    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    /*constexpr double x1 = -PI;
    constexpr double x2 = PI;
    constexpr double y1 = -PI;
    constexpr double y2 = PI;

    constexpr double Ho = 0.2;
    constexpr double zo = 0.025;
    Utilities::ignore(Ho);

    constexpr double w   = 1;
    constexpr double tau = 0;

    double ic_ze = 2 * zo * cos(w * (x - x1)) * cos(w * (y - y1)) *
                   cos(w * (t + tau)) / (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double ic_qx = zo * sin(w * (x - x1)) * cos(w * (y - y1)) * sin(w * (t + tau)) /
                   (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double ic_qy = zo * cos(w * (x - x1)) * sin(w * (y - y1)) * sin(w * (t + tau)) /
                   (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));*/

    double ic_ze = exp(sin(3 * x) * sin(3 * y) - sin(3 * t));
    double ic_qx = cos(x - 4 * t);
    double ic_qy = sin(y + 4 * t);

    StatVector<double, SWE::n_variables> ic_q{ic_ze, ic_qx, ic_qy};

    return ic_q;
}
}