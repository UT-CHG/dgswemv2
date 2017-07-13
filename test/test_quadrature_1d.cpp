#include <cmath>
#include <iostream>
#include <limits>

#include "geometry/quadrature/gauss_legendre.hpp"

namespace Geometry {
  namespace Basis {
    std::vector<double> jacobi_polynomial(uint n, int a, int b, const std::vector<double>& x);
  }
}

using Geometry::Basis::jacobi_polynomial;

int main() {

  bool any_error = false;

  for ( uint p = 0; p < 20; ++p ) {

    bool error_found = false;

    Geometry::Quadrature::GaussLegendre quadr(p);

    //set gps
    std::vector<double> gps=quadr.get_gp();

    for ( uint q =0; q <= p; ++q ) {
      std::vector<double> legendre_evals = jacobi_polynomial(q, 0, 0, gps);

      double eval = quadr.integrate(legendre_evals);

      if ( q == 0 ) {
        if ( std::abs(eval-2) > 100*std::numeric_limits<double>::epsilon() ) {
          error_found = true;
        }
      } else {
        if ( std::abs(eval) > 100*std::numeric_limits<double>::epsilon() ) {
          error_found = true;
        }
      }
    }

    if ( error_found ) {
      std::cerr << "Error found in GaussLegendre at " << std::to_string(p) << std::endl;
      any_error = true;
    }
  }

  if ( any_error ) {
    return 1;
  }

 return 0;

}
