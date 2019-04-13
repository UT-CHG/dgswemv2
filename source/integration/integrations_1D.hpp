#ifndef INTEGRATIONS_1D_HPP
#define INTEGRATIONS_1D_HPP

#include "general_definitions.hpp"

namespace Integration {
/**
 * Gauss Legendre quadrature rule for the unit interval.
 * GaussLegendre_1D describes the Gauss Legendre quadrature rules for the unit interval.
 * @note The strength of the rule is limited to polynomial degree 65.
 */
class GaussLegendre_1D : Integration<1> {
  public:
    std::pair<DynVector<double>, AlignedVector<Point<1>>> GetRule(const uint p);

    uint GetNumGP(const uint p);

  private:
    /**
     * Get the data for the Gauss Legendre quadrature rules.
     *
     * @param number of gauss points desired in the rule.
     * @return Pair of the weights and point on the unit interval
     */
    std::pair<std::vector<double>, AlignedVector<Point<1>>> GPData(const uint number_gp);
};
}

#endif