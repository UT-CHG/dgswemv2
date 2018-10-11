#ifndef BASES_1D_HPP
#define BASES_1D_HPP

#include "general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
class Legendre_1D : Basis<1> {
  public:
    DynMatrix<double> GetPhi(const uint p, const std::vector<Point<1>>& points);
    std::array<DynMatrix<double>, 1> GetDPhi(const uint p, const std::vector<Point<1>>& points);

    DynMatrix<double> GetMinv(const uint p);

    DynMatrix<double> ProjectBasisToLinear(const DynMatrix<double>& u);
    DynMatrix<double> ProjectLinearToBasis(const uint ndof, const DynMatrix<double>& u_lin);
};
}

#endif