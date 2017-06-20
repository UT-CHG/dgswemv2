#ifndef BASIS_POLYNOMIALS_H
#define BASIS_POLYNOMIALS_H

#include "../general_definitions.h"

namespace Basis {
	std::vector<double> jacobi_polynomial(int n, int a, int b, const std::vector<double>& x);
	std::vector<double> jacobi_polynomial_derivative(int n, int a, int b, const std::vector<double>& x);
}

#endif