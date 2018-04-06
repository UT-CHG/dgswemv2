#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    double bath;

    if (x < 255000.0) {
        bath = (x - 225000.0) / 15000.0;
    } else {
        bath = 50.0 - 50.0 * exp(-(x - 225000.0) / 750000.0);
    }

    return bath;  // bathymetry for storm
}

#endif