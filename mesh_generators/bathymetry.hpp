#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    double bath = 0.0;

    if (std::hypot(x - 2.0, y) >= 2.0) {
        bath = 0.3;
    } else if (std::hypot(x - 2.0, y) < 2.0) {
        bath = 0.3 - 1.0 / 8.0 * (1 + std::cos(PI / 2.0 * std::hypot(x - 2.0, y)));
    }

    return bath;
}
#endif
