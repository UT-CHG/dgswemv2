#ifndef BASES_2D_HPP
#define BASES_2D_HPP

#include "../general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
class Dubiner_2D : Basis<2> {
  public:
    Array2D<double> GetPhi(const uint p, const std::vector<Point<2>>& points);
    Array3D<double> GetDPhi(const uint p, const std::vector<Point<2>>& points);
    std::pair<bool, Array2D<double>> GetMinv(const uint p);

  private:
    std::vector<double> ComputePhi(const uint p,
                                   const uint q,
                                   const std::vector<double>& n1,
                                   const std::vector<double>& n2);
    Array2D<double> ComputeDPhi(const uint p,
                                const uint q,
                                const std::vector<double>& n1,
                                const std::vector<double>& n2);

    std::vector<double> ComputeSingularDPhi(const uint p, const uint q);
    std::vector<double> ComputeSingularDPhiDZ1(const uint q);
    std::vector<double> ComputeSingularDPhiDZ2(const uint q);
};
}

#endif