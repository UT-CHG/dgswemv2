#ifndef BASIS_POLYNOMIALS_HPP
#define BASIS_POLYNOMIALS_HPP

#include "../general_definitions.hpp"

namespace Basis {
std::vector<double> jacobi_polynomial(uint n, uint a, uint b, const std::vector<double>& x);
std::vector<double> jacobi_polynomial_derivative(uint n, uint a, uint b, const std::vector<double>& x);
}

#endif