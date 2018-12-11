#include "problem/SWE/problem_function_files/fwd.hpp"

namespace SWE {
StatVector<double, SWE::n_variables> ic_q(const double t, const Point<2>& pt) {
    double ic_ze = 0.0;
    double ic_qx = 0.0;
    double ic_qy = 0.0;

    StatVector<double, SWE::n_variables> ic_q{ic_ze, ic_qx, ic_qy};

    return ic_q;
}
}