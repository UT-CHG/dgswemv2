#include <cmath>
#include <iostream>
#include <limits>
#include <memory>


#include "geometry/basis/dubiner.hpp"
#include "geometry/quadrature/dunavant.hpp"

typedef std::array<double,2> Point;

bool check_quadrature( int p,
                       const Geometry::Quadrature::Base<Point>* quadr,
                       Geometry::Basis::Dubiner& basis)
{
  auto gps = quadr->get_gp();
  bool error_found = false;

  for ( uint dof =0; dof < (p+1)*(p+2)/2; ++dof ) {
    auto dof_evals = basis.get(dof,gps);

    double eval = quadr->integrate(dof_evals);

    if ( dof == 0 ) {
      if ( std::abs(eval-2) > 100*std::numeric_limits<double>::epsilon() ) {
        error_found = true;
      }
    } else {
      if ( std::abs(eval) > 100*std::numeric_limits<double>::epsilon() ) {
        error_found = true;
      }
    }
  }

  return error_found;
}

int main() {

  Geometry::Basis::Dubiner basis(20);
  bool any_error = false;

  for ( uint p = 1; p < 20; ++p ) {

    std::unique_ptr<Geometry::Quadrature::Base<Point> > quadr( new Geometry::Quadrature::Dunavant(p) );
    bool error = check_quadrature(p, quadr.get(),
                                  basis);

    if ( error ) {
      any_error = true;
      std::cerr << "Error found in Dubiner at " << std::to_string(p) << std::endl;
    }
  }

  if ( any_error ) {
    return 1;
  }

  return 0;
}
