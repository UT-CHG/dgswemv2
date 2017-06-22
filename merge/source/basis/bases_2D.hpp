#ifndef BASES_2D_HPP
#define BASES_2D_HPP

#include "../general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
	class Dubiner_2D : Basis<2> {
	public:
		Array2D<double> get_phi(uint, const std::vector<Point<2>>&);
		Array3D<double> get_dphi(uint, const std::vector<Point<2>>&);
		std::pair<bool, Array2D<double>> get_m_inv(uint);
		void basis_test(uint, const Array2D<double>&, const std::pair<std::vector<double>, std::vector<Point<2>>>&);
	private:
		std::vector<double> dubiner_2d_phi(uint, uint, const std::vector<double>&, const std::vector<double>&);
		Array2D<double> dubiner_2d_dphi(uint, uint, const std::vector<double>&, const std::vector<double>&);
	};
}

#endif