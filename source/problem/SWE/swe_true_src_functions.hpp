#ifndef SWE_TRUE_SRC_FUNCTIONS_HPP
#define SWE_TRUE_SRC_FUNCTIONS_HPP

namespace SWE {
inline double true_ze(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    return 2 * zo * cos(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) * cos(w * (t + tau)) /
           (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));
}

inline double true_qx(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    return zo * sin(w * (pt[GlobalCoord::x] - x1)) * cos(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
           (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));
}

inline double true_qy(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    return zo * cos(w * (pt[GlobalCoord::x] - x1)) * sin(w * (pt[GlobalCoord::y] - y1)) * sin(w * (t + tau)) /
           (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));
}

inline double source_ze(const double t, const Point<2>& pt) { return 0; }

inline double source_qx(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    return w * zo * cos((t + tau) * w) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))) *
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
               pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) *
                            (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))),
                   2) +
           (2 * w * pow(zo, 3) * cos((t + tau) * w) * pow(cos(w * (x - x1)), 2) * cos(w * (y - y1)) *
            pow((1. / cos(w * (-x1 + x2))), 3) * pow((1. / cos(w * (-y1 + y2))), 3) * pow(sin((t + tau) * w), 2) *
            sin(w * (x - x1)) * pow(sin(w * (y - y1)), 2)) /
               pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) *
                            (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))),
                   2) -
           (w * pow(zo, 2) * cos(w * (x - x1)) * pow((1. / cos(w * (-x1 + x2))), 2) *
            pow((1. / cos(w * (-y1 + y2))), 2) * pow(sin((t + tau) * w), 2) * sin(w * (x - x1)) *
            pow(sin(w * (y - y1)), 2)) /
               (Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                         (1. / cos(w * (-y1 + y2))));
}

inline double source_qy(const double t, const Point<2>& pt) {
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;

    constexpr double Ho = 2.;
    constexpr double zo = 0.25;

    constexpr double w = 2 * PI / 43200.;
    constexpr double tau = 0;

    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    return w * zo * cos((t + tau) * w) * cos(w * (x - x1)) * (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))) *
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
               pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) *
                            (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))),
                   2) -
           (w * pow(zo, 2) * cos(w * (y - y1)) * pow((1. / cos(w * (-x1 + x2))), 2) *
            pow((1. / cos(w * (-y1 + y2))), 2) * pow(sin((t + tau) * w), 2) * pow(sin(w * (x - x1)), 2) *
            sin(w * (y - y1))) /
               (Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) * (1. / cos(w * (-x1 + x2))) *
                         (1. / cos(w * (-y1 + y2)))) +
           (2 * w * pow(zo, 3) * cos((t + tau) * w) * pow(cos(w * (x - x1)), 3) * pow((1. / cos(w * (-x1 + x2))), 3) *
            pow((1. / cos(w * (-y1 + y2))), 3) * pow(sin((t + tau) * w), 2) * pow(sin(w * (y - y1)), 3)) /
               pow(Ho + 2 * zo * cos((t + tau) * w) * cos(w * (x - x1)) * cos(w * (y - y1)) *
                            (1. / cos(w * (-x1 + x2))) * (1. / cos(w * (-y1 + y2))),
                   2);
}
}

#endif