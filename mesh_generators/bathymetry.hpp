#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    double bath = 0;

    if (x <= -80) {
        bath = 60.0 * (x + 95.0);
    } else {
        bath = 900.0;
    }

    return bath;  // bathymetry for storm
}

#endif