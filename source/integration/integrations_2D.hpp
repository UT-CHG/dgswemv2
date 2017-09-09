#ifndef INTEGRATIONS_2D_HPP
#define INTEGRATIONS_2D_HPP

#include "../general_definitions.hpp"

namespace Integration {
/**
 * Quadrature rule for the master triangle.
 * Dunavant_2D describes a numerical quadrature rule for a triangle with vertices at
 * (-1,1), (1,-1), and (-1,1)
 */
class Dunavant_2D : Integration<2> {
  public:
    /**
     * Returns set of weights and quadrature points.
     * GetRule returns a vector of quadrature points and weights over
     * a triangle with vertices at (-1,-1), (1,-1), and (-1,1).
     * @param p Polynomial order for which the rule should return exact results.
     * \note `p` cannot exceed 20.
     */
    std::pair<std::vector<double>, std::vector<Point<2>>> GetRule(const uint p);
    /**
     * Returns the number of Gauss points required for a rule of strength p.
     * GetNumGP returns the number of quadrature points for a rule of strength p
     * without having to construct the rule.
     */
    uint GetNumGP(const uint p)

  private:
    std::vector<uint> PermutationData(const uint p);
    std::pair<std::vector<double>, std::vector<Point<3>>> GPData(const uint p);
};
}

#endif
