#ifndef CLASS_MASTER_ELEMENT_HPP
#define CLASS_MASTER_ELEMENT_HPP

#include "general_definitions.hpp"

namespace Master {
/**
 * Triangular master element.
 * Triangular master elements contains the information necessary for evaluating computations on
 * on the master triangle.
 */
template <typename BasisType, typename IntegrationType>
class Triangle : public Master<2> {
  public:
    /**
     * The basis used over the element.
     */
    BasisType basis;

    /**
     * Integration rule used for evaluating integrals over the element.
     */
    IntegrationType integration;

  public:
    /**
     * Default constructor
     */
    Triangle() = default;

    /**
     * Construct a master triangle with polynomial order p.
     *
     * @param p Polynomial order
     */
    Triangle(const uint p);

    /**
     * Transform coordinates on a boundary to master element coordinates.
     * The function accepts a boundary ID and a vector of points located on that boundary
     * and returns their coordinates on the master triangle. Based on boundary ID the endpoints
     * are mapped to the following coordinates:
     *
     * Boundary ID | -1 Mapped To: | 1 Mapped To:
     * ------------|---------------|-------------
     *      0      |   ( 1,-1)     |   (-1, 1)
     *      1      |   (-1, 1)     |   (-1,-1)
     *      2      |   (-1,-1)     |   ( 1,-1)
     *
     * @param bound_id The ID of the boundary
     * @param The points on the boundary
     */
    DynVector<Point<2>> BoundaryToMasterCoordinates(const uint bound_id, const DynVector<Point<1>>& z_boundary);

    template <typename InputArrayType>
    decltype(auto) ComputeLinearUbaryctr(const InputArrayType& u_lin);
    template <typename InputArrayType>
    decltype(auto) ComputeLinearUmidpts(const InputArrayType& u_lin);
    template <typename InputArrayType>
    decltype(auto) ComputeLinearUvrtx(const InputArrayType& u_lin);

  private:
    DynVector<Point<2>> VTKPostCell();
    DynVector<Point<2>> VTKPostPoint();
};
}

#include "elements_2D/master_triangle.tpp"

#endif
