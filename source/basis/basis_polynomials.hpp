#ifndef BASIS_POLYNOMIALS_HPP
#define BASIS_POLYNOMIALS_HPP

#include "../general_definitions.hpp"

namespace Basis {
std::vector<double> jacobi_polynomial(const uint n, const uint a, const uint b, const std::vector<double>& x);
std::vector<double> jacobi_polynomial_derivative(const uint n,
                                                 const uint a,
                                                 const uint b,
                                                 const std::vector<double>& x);
}

#endif