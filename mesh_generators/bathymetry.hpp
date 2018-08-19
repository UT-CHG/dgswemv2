#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    double bath = 0.0;

    if (x < 20.0) {
        bath = 0.7;
    } else if (x >= 20.0) {
        bath = 0.7 - 0.02 * (x - 20.0);
    }

    return bath;
}
#endif
