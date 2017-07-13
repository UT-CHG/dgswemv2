#ifndef BASES_2D_HPP
#define BASES_2D_HPP

#include "../general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
	class Dubiner_2D : Basis<2> {
	public:
		Array2D<double> GetPhi(uint, const std::vector<Point<2>>&);
		Array3D<double> GetDPhi(uint, const std::vector<Point<2>>&);
		std::pair<bool, Array2D<double>> GetMinv(uint);
	private:
		std::vector<double> ComputePhi(uint, uint, const std::vector<double>&, const std::vector<double>&);
		Array2D<double> ComputeDPhi(uint, uint, const std::vector<double>&, const std::vector<double>&);

		std::vector<double> ComputeSingularDPhi(uint, uint);
		double ComputeSingularDPhiDZ1(uint);
		double ComputeSingularDPhiDZ2(uint);
	};
}

#endif