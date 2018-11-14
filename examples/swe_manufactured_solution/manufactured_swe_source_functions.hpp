#ifndef SWE_SOURCE_FUNCTIONS_HPP
#define SWE_SOURCE_FUNCTIONS_HPP

namespace SWE {
inline StatVector<double, SWE::n_variables> source_q(const double t, const Point<2>& pt) {
    /*constexpr double x1 = -PI;
    constexpr double x2 = PI;
    constexpr double y1 = -PI;
    constexpr double y2 = PI;

    constexpr double Ho = 0.2;
    constexpr double zo = 0.025;

    constexpr double w   = 1;
    constexpr double tau = 0;

    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    double source_ze = 0.0;

    double source_qx =
        w * zo * cos((t + tau) * w) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))) *
            sin(w * (x - x1)) -
        2 * SWE::Global::g * Ho * w * zo * cos((t + tau) * w) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
            (1. / cos(w * (-y1 + y2))) * sin(w * (x - x1)) -
        4. * SWE::Global::g * w * pow(zo, 2) * pow(cos((t + tau) * w), 2) * cos(w * (x - x1)) *
            pow(cos(w * (y - y1)), 2) * pow((1. / cos(w * (-x1 + x2))), 2) * pow((1. / cos(w * (-y1 + y2))), 2) *
            sin(w * (x - x1)) +
        (3 * w * pow(zo, 2) * cos(w * (x - x1)) * pow(cos(w * (y - y1)), 2) * pow((1. / cos(w * (-x1 + x2))), 2) *
         pow((1. / cos(w * (-y1 + y2))), 2) * pow(sin((t + tau) * w), 2) * sin(w * (x - x1))) /
            (Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                      (1. / cos(w * (-y1 + y2)))) +
        (2 * w * pow(zo, 3) * cos((t + tau) * w) * pow(cos(w * (y - y1)), 3) * pow((1. / cos(w * (-x1 + x2))), 3) *
         pow((1. / cos(w * (-y1 + y2))), 3) * pow(sin((t + tau) * w), 2) * pow(sin(w * (x - x1)), 3)) /
            pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                         (1. / cos(w * (-y1 + y2))),
                2) +
        (2 * w * pow(zo, 3) * cos((t + tau) * w) * pow(cos(w * (x - x1)), 2) * cos(w * (y - y1)) *
         pow((1. / cos(w * (-x1 + x2))), 3) * pow((1. / cos(w * (-y1 + y2))), 3) * pow(sin((t + tau) * w), 2) *
         sin(w * (x - x1)) * pow(sin(w * (y - y1)), 2)) /
            pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                         (1. / cos(w * (-y1 + y2))),
                2) -
        (w * pow(zo, 2) * cos(w * (x - x1)) * pow((1. / cos(w * (-x1 + x2))), 2) * pow((1. / cos(w * (-y1 + y2))), 2) *
         pow(sin((t + tau) * w), 2) * sin(w * (x - x1)) * pow(sin(w * (y - y1)), 2)) /
            (Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                      (1. / cos(w * (-y1 + y2))));

    double source_qy =
        w * zo * cos((t + tau) * w) * cos(w * (x - x1)) * (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))) *
            sin(w * (y - y1)) -
        2 * SWE::Global::g * Ho * w * zo * cos((t + tau) * w) * cos(w * (x - x1)) * (1. / cos(w * (-x1 + x2))) *
            (1. / cos(w * (-y1 + y2))) * sin(w * (y - y1)) -
        4. * SWE::Global::g * w * pow(zo, 2) * pow(cos((t + tau) * w), 2) * pow(cos(w * (x - x1)), 2) *
            cos(w * (y - y1)) * pow((1. / cos(w * (-x1 + x2))), 2) * pow((1. / cos(w * (-y1 + y2))), 2) *
            sin(w * (y - y1)) +
        (3 * w * pow(zo, 2) * pow(cos(w * (x - x1)), 2) * cos(w * (y - y1)) * pow((1. / cos(w * (-x1 + x2))), 2) *
         pow((1. / cos(w * (-y1 + y2))), 2) * pow(sin((t + tau) * w), 2) * sin(w * (y - y1))) /
            (Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                      (1. / cos(w * (-y1 + y2)))) +
        (2 * w * pow(zo, 3) * cos((t + tau) * w) * cos(w * (x - x1)) * pow(cos(w * (y - y1)), 2) *
         pow((1. / cos(w * (-x1 + x2))), 3) * pow((1. / cos(w * (-y1 + y2))), 3) * pow(sin((t + tau) * w), 2) *
         pow(sin(w * (x - x1)), 2) * sin(w * (y - y1))) /
            pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                         (1. / cos(w * (-y1 + y2))),
                2) -
        (w * pow(zo, 2) * cos(w * (y - y1)) * pow((1. / cos(w * (-x1 + x2))), 2) * pow((1. / cos(w * (-y1 + y2))), 2) *
         pow(sin((t + tau) * w), 2) * pow(sin(w * (x - x1)), 2) * sin(w * (y - y1))) /
            (Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                      (1. / cos(w * (-y1 + y2)))) +
        (2 * w * pow(zo, 3) * cos((t + tau) * w) * pow(cos(w * (x - x1)), 3) * pow((1. / cos(w * (-x1 + x2))), 3) *
         pow((1. / cos(w * (-y1 + y2))), 3) * pow(sin((t + tau) * w), 2) * pow(sin(w * (y - y1)), 3)) /
            pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                         (1. / cos(w * (-y1 + y2))),
                2);*/

    double source_ze = -3 * exp(-sin(3 * t) + sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])) * cos(3 * t) +
                       cos(4 * t + pt[GlobalCoord::y]) + sin(4 * t - pt[GlobalCoord::x]);

    double source_qx =
        (-3 * exp(3 * sin(3 * t) + sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])) *
             pow(cos(4 * t - pt[GlobalCoord::x]), 2) * cos(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]) +
         (2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) *
             (exp(3 * sin(3 * t)) * sin(8 * t - 2 * pt[GlobalCoord::x]) -
              (2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) *
                  (4 * exp(2 * sin(3 * t)) * sin(4 * t - pt[GlobalCoord::x]) -
                   3 * exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])) *
                       (2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) *
                       cos(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) +
         exp(3 * sin(3 * t)) * cos(4 * t - pt[GlobalCoord::x]) *
             ((2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) *
                  cos(4 * t + pt[GlobalCoord::y]) -
              3 * exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])) * cos(3 * pt[GlobalCoord::y]) *
                  sin(3 * pt[GlobalCoord::x]) * sin(4 * t + pt[GlobalCoord::y]))) /
        (exp(2 * sin(3 * t)) *
         pow(2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])), 2));

    double source_qy =
        (4 * exp(2 * sin(3 * t)) *
             pow(2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])), 2) *
             cos(4 * t + pt[GlobalCoord::y]) +
         3 * exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])) * cos(3 * pt[GlobalCoord::y]) *
             sin(3 * pt[GlobalCoord::x]) *
             (pow(2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])), 3) -
              exp(3 * sin(3 * t)) * pow(sin(4 * t + pt[GlobalCoord::y]), 2)) +
         exp(3 * sin(3 * t)) *
             ((2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) *
                  sin(4 * t - pt[GlobalCoord::x]) * sin(4 * t + pt[GlobalCoord::y]) -
              3 * exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])) * cos(4 * t - pt[GlobalCoord::x]) *
                  cos(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]) * sin(4 * t + pt[GlobalCoord::y]) +
              (2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y]))) *
                  sin(2 * (4 * t + pt[GlobalCoord::y])))) /
        (exp(2 * sin(3 * t)) *
         pow(2 * exp(sin(3 * t)) + exp(sin(3 * pt[GlobalCoord::x]) * sin(3 * pt[GlobalCoord::y])), 2));

    StatVector<double, SWE::n_variables> source_q{source_ze, source_qx, source_qy};

    return source_q;
}
}

#endif
