#ifndef BASES_1D_HPP
#define BASES_1D_HPP

#include "general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
class Legendre_1D : Basis<1> {
  public:
    DynMatrix<double> GetPhi(const uint p, const AlignedVector<Point<1>>& points) override;
    std::array<DynMatrix<double>, 1> GetDPhi(const uint p, const AlignedVector<Point<1>>& points) override;

    DynMatrix<double> GetMinv(const uint p) override;

    DynMatrix<double> GetBasisLinearT(const uint p) override;
    DynMatrix<double> GetLinearBasisT(const uint p) override;
};
}

#endif