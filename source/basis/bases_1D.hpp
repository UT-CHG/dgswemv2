#ifndef BASES_1D_HPP
#define BASES_1D_HPP

#include "general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
class Legendre_1D : Basis<1> {
  public:
    DynMatrix<double> GetPhi(const uint p, const DynVector<Point<1>>& points);
    StatVector<DynMatrix<double>, 1> GetDPhi(const uint p, const DynVector<Point<1>>& points);

    std::pair<bool, DynMatrix<double>> GetMinv(const uint p);

    template <typename T>
    void ProjectBasisToLinear(const std::vector<T>& u, std::vector<T>& u_lin);
    template <typename T>
    void ProjectLinearToBasis(const std::vector<T>& u_lin, std::vector<T>& u);
};
}

#endif