#ifndef BASES_2D_HPP
#define BASES_2D_HPP

#include "../general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
/**
 * Dubiner basis class.
 * This class implements a basis over the reference triangle given via coordinates
 * (-1,-1), (1,-1), and (-1,1).
 */
class Dubiner_2D : Basis<2> {
  public:
    Array2D<double> GetPhi(const uint p, const std::vector<Point<2>>& points);
    Array3D<double> GetDPhi(const uint p, const std::vector<Point<2>>& points);

    std::pair<bool, Array2D<double>> GetMinv(const uint p);

    void ProjectBasisToLinear(const std::vector<double>& u, std::vector<double>& u_lin);
    void ProjectLinearToBasis(const std::vector<double>& u_lin, std::vector<double>& u);

  private:
    /**
     * Compute the (p,q)-th dubiner polynomial at points (n1,n2).
     *
     * @param p
     * @param q
     * @param n1 set of x-coordinates to evaluate the Dubiner polynomials at.
     * @param n2 set of y-corrdinates to evaluate the Dubiner polynomials at.
     * @return evaluations of the Dubiner polynomials at points (n1,n2)
     * @warning It is required that the size of n1 is equal to the size of n2.
     */
    std::vector<double> ComputePhi(const uint                 p,
                                   const uint                 q,
                                   const std::vector<double>& n1,
                                   const std::vector<double>& n2);

    /**
     * Compute the gradient of the (p,q)-th dubiner polynomial at points (n1,n2).
     *
     * @param p
     * @param q
     * @param n1 set of x-coordinates to evaluate the Dubiner polynomials at.
     * @param n2 set of y-corrdinates to evaluate the Dubiner polynomials at.
     * @return evaluations of the gradient at points (n1,n2)
     * @warning It is required that the size of n1 is equal to the size of n2.
     */
    Array2D<double> ComputeDPhi(const uint                 p,
                                const uint                 q,
                                const std::vector<double>& n1,
                                const std::vector<double>& n2);

    /**
     * Get of the gradient of the Dubiner polynomial at the singularity.
     * Evaluate of the gradient of the (p,q)-th Dubiner polynomial at the
     * singularity at (-1,1).
     *
     * @param p
     * @param q
     * @return the gradient of the (p,q)-th Dubiner polynomial
     */
    std::vector<double> ComputeSingularDPhi(const uint p, const uint q);

    /**
     * Get the derivative of Dubiner polynomial in the 1st coordinate at the singularity.
     * Evaluate the derivative in the 1st coordinate of the (p,q)-th Dubiner polynomial at the
     * singularity at (-1,1).
     *
     * @param p
     * @param q
     * @return the derivative in the 1st coordinate of the (p,q)-th Dubiner polynomial
     */
    std::vector<double> ComputeSingularDPhiDZ1(const uint q);

    /**
     * Get the derivative of Dubiner polynomial in the 2nd coordinate at the singularity.
     * Evaluate the derivative in the 2nd coordinate of the (p,q)-th Dubiner polynomial at the
     * singularity at (-1,1).
     *
     * @param p
     * @param q
     * @return the derivative in the 2nd coordinate of the (p,q)-th Dubiner polynomial
     */
    std::vector<double> ComputeSingularDPhiDZ2(const uint q);
};
}

#endif