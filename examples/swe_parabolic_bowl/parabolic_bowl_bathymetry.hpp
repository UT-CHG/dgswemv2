#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    constexpr double alpha = 1.e-7;

    return -alpha * (x*x + y*y);
}
#endif
