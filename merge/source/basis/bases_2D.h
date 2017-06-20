#ifndef BASES_2D_H
#define BASES_2D_H

#include "../general_definitions.h"
#include "basis_polynomials.h"

namespace Basis {
	class Dubiner_2D : Basis<2> {
	public:
		Array2D<double> get_phi(int, const std::vector<Point<2>>&);
		Array3D<double> get_dphi(int, const std::vector<Point<2>>&);
		std::pair<bool, Array2D<double>> get_m_inv(int);
		void basis_test(int, const Array2D<double>&, const std::pair<std::vector<double>, std::vector<Point<2>>>&);
	private:
		std::vector<double> dubiner_2d_phi(int, int, const std::vector<double>&, const std::vector<double>&);
		Array2D<double> dubiner_2d_dphi(int, int, const std::vector<double>&, const std::vector<double>&);
	};
}

#endif