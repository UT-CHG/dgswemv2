#ifndef INTEGRATIONS_2D_HPP
#define INTEGRATIONS_2D_HPP

#include "../general_definitions.hpp"

namespace Integration {
/**
 * Quadrature rule for the master triangle.
 * Dunavant_2D describes a numerical quadrature rule for a triangle with vertices at
 * (-1,1), (1,-1), and (-1,1)
 * @note The strength of the rule is limited to polynomial degree 20.
 */
class Dunavant_2D : Integration<2> {
  public:
    // Implementations of inherited functions
    std::pair<std::vector<double>, std::vector<Point<2>>> GetRule(const uint p);
    uint GetNumGP(const uint p);

  private:
    /**
     * Get the size of the permutations of gauss points under the symmetries of the triangle.
     * The quadrature rule can be expressed compactly by modding out by reflections
     * and rotations of the triangle. The permutation size of a specific gauss point
     *  corresponds to the size of its orbit under D<sub>3</sub>.
     *
     * @param p Polynomial order for which the rule should return exact results.
     * @return Vector of permutation sizes.
     */
    std::vector<uint> PermutationData(const uint p);

    /**
     * Get the compact representation of the quadrature rule.
     * The Dunavant basis is independent of the orientation of the triangle as such
     * it is possible to sparsely represent the coordinates by modding out symmetric points.
     * GPData obtains this compact representation of the quadrature rule
     *
     * @param p Polynomial order for which the rule should return exact results.
     * @return Pair of weights and coordinates in Barycentric coordinates.
     */
    std::pair<std::vector<double>, std::vector<Point<3>>> GPData(const uint p);
};

class GaussLegendre_2D : Integration<2> {
  public:
    std::pair<std::vector<double>, std::vector<Point<2>>> GetRule(const uint p);
    uint GetNumGP(const uint p);

  private:
  /**trying to figure out what the hell I'm doing right now
   * 
   * 
   * 
   * 
   *
   */

};

}

#endif
