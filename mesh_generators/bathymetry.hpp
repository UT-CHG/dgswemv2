#ifndef BATHYMETRY_HPP
#define BATHYMETRY_HPP

double bathymetry_function(double x, double y) {
    double bath = 3.0 - (x / 10000.0) * (x / 10000.0) - (y / 10000.0) * (y / 10000.0);

    return bath;
}

#endif