#ifndef BASES_2D_HPP
#define BASES_2D_HPP

#include "general_definitions.hpp"
#include "basis_polynomials.hpp"

namespace Basis {
/**
 * Dubiner basis class.
 * This class implements a basis over the reference triangle given via coordinates
 * (-1,-1), (1,-1), and (-1,1).
 */
class Dubiner_2D : Basis<2> {
  public:
    DynMatrix<double> GetPhi(const uint p, const std::vector<Point<2>>& points);
    std::array<DynMatrix<double>, 2> GetDPhi(const uint p, const std::vector<Point<2>>& points);

    DynMatrix<double> GetMinv(const uint p);

    DynMatrix<double> ProjectBasisToLinear(const DynMatrix<double>& u);
    DynMatrix<double> ProjectLinearToBasis(const uint ndof, const DynMatrix<double>& u_lin);

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
    DynVector<double> ComputePhi(const uint p,
                                 const uint q,
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
    DynVector<double> ComputeDPhiDZ1(const uint p,
                                     const uint q,
                                     const std::vector<double>& n1,
                                     const std::vector<double>& n2);

    DynVector<double> ComputeDPhiDZ2(const uint p,
                                     const uint q,
                                     const std::vector<double>& n1,
                                     const std::vector<double>& n2);

    /**
     * Get the derivative of Dubiner polynomial in the 1st coordinate at the singularity.
     * Evaluate the derivative in the 1st coordinate of the (p,q)-th Dubiner polynomial at the
     * singularity at (-1,1).
     *
     * @param p
     * @param q
     * @return the derivative in the 1st coordinate of the (p,q)-th Dubiner polynomial
     */
    std::array<double, 2> ComputeSingularDPhiDZ1(const uint q);

    /**
     * Get the derivative of Dubiner polynomial in the 2nd coordinate at the singularity.
     * Evaluate the derivative in the 2nd coordinate of the (p,q)-th Dubiner polynomial at the
     * singularity at (-1,1).
     *
     * @param p
     * @param q
     * @return the derivative in the 2nd coordinate of the (p,q)-th Dubiner polynomial
     */
    std::array<double, 2> ComputeSingularDPhiDZ2(const uint q);
};
}

#endif