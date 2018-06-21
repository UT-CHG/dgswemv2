#ifndef BASES_1D_HPP
#define BASES_1D_HPP

#include "general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
class Legendre_1D : Basis<1> {
  public:
    Array2D<double> GetPhi(const uint p, const std::vector<Point<1>>& points);
    Array3D<double> GetDPhi(const uint p, const std::vector<Point<1>>& points);

    std::pair<bool, Array2D<double>> GetMinv(const uint p);

    void ProjectBasisToLinear(const std::vector<double>& u, std::vector<double>& u_lin);
    void ProjectLinearToBasis(const std::vector<double>& u_lin, std::vector<double>& u);
};
}

#endif