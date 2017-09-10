#ifndef BASIS_POLYNOMIALS_HPP
#define BASIS_POLYNOMIALS_HPP

#include "../general_definitions.hpp"

namespace Basis {
/**
 * Evaluate the Jacobi polynomial P<sub>n</sub><sup>(a,b)</sup> at points x
 *
 * @param n
 * @param a
 * @param b
 * @param x The points at which the Jacobi polynomial is evaluated.
 * @return Jacobi polynomial evaluation at the points x
 */
std::vector<double> jacobi_polynomial(const uint n, const uint a, const uint b, const std::vector<double>& x);

/**
 * Evaluate the derivative of the Jacobi polynomial P<sub>n</sub><sup>(a,b)</sup> at points x
 *
 * @param n
 * @param a
 * @param b
 * @param x The points at which the Jacobi polynomial is evaluated.
 * @return Jacobi polynomial's derivate evaluation at the points x
 */
std::vector<double> jacobi_polynomial_derivative(const uint n,
                                                 const uint a,
                                                 const uint b,
                                                 const std::vector<double>& x);
}

#endif