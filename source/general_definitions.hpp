#ifndef GENERAL_DEFINITIONS
#define GENERAL_DEFINITIONS

#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <assert.h>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>

typedef unsigned int uint;

typedef unsigned char uchar;

template <uint dimension>
using Point = std::array<double, dimension>;

template <class type>
using Array2D = std::vector<std::vector<type>>;

template <class type>
using Array3D = std::vector<std::vector<std::vector<type>>>;

template <class type>
using Array4D = std::vector<std::vector<std::vector<std::vector<type>>>>;

namespace Basis {
/**
 * Base class for the basis types over the element.
 *
 * @tparam dimension over the surface to be integrated
 */
template <uint dimension>
class Basis {
  public:
    /**
     * Evaluate the basis function at points.
     * Evaluate all basis functions evaluated at points up to polynomial order p.
     * e.g. in the one dimensional case, p+1 functions will be evaluated, and in the two
     * dimensional case (p+1)(p+2)/2 functions will be evaluated.
     *
     * @param p polynomial order
     * @param points vector of points at which the basis functions are evaluated
     * @return 2-dimensional array indexed as [dof][point]
     */
    virtual Array2D<double> GetPhi(const uint p, const std::vector<Point<dimension>>& points) = 0;
    /**
     * Evaluate the gradients of the basis function at points.
     * Evaluate gradients of all basis functions evaluated at points up to polynomial order p.
     * e.g. in the one dimensional case, p+1 functions will be evaluated, and in the two
     * dimensional case (p+1)(p+2)/2 functions will be evaluated.
     *
     * @param p polynomial order
     * @param points vector of points at which the basis functions are evaluated
     * @return 3-dimensional array indexed as [dof][point][coordinate] where coordinate
     *         corresponds to the dimension in which the derivative is taken.
     */
    virtual Array3D<double> GetDPhi(const uint p, const std::vector<Point<dimension>>& points) = 0;

    /**
     * Obtain the inverted mass matrix of basis functions of polynomial order p.
     *
     * @param p polynomial order
     * @param return a pair with a boolean, which states whether the basis is orthogonal,
     *        and a 2-dimensional array corresponding to the mass matrix over the master element
     */
    virtual std::pair<bool, Array2D<double>> GetMinv(const uint p) = 0;
};
}

namespace Integration {
/**
 * Base class for Integration types over the element.
 *
 * @tparam dimension over the surface to be integrated
 */
template <uint dimension>
class Integration {
  public:
    /**
     * Returns the number of Gauss points required for a rule of strength p
     * GetNumGP returns the number of quadrature points for a rule of strength p without
     * without having to construct the rule.
     *
     * @param p Polynomial order for which the rule should return exact results.
     * @return Pair of weights and quadrature points on the master triangle.
     */
    virtual uint GetNumGP(const uint p) = 0;

    /**
     * Returns a vector of weights and quadrature points.
     *
     * @param p Polynomial order for which the rule should return exact results.
     * @return Number of Gauss points
     */
    virtual std::pair<std::vector<double>, std::vector<Point<dimension>>> GetRule(const uint p) = 0;
};
}

namespace Master {
template <uint dimension>
class Master {
  public:
    uint p;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule;

    Array2D<double> phi_gp;
    Array3D<double> dphi_gp;

    Array2D<double> int_fact_phi;
    Array3D<double> int_fact_dphi;

    std::pair<bool, Array2D<double>> m_inv;

    Array2D<double> phi_postprocessor_cell;
    Array2D<double> phi_postprocessor_point;

  public:
    Master(const uint p) : p(p) {}

    virtual std::vector<Point<dimension>> BoundaryToMasterCoordinates(
        const uint bound_id,
        const std::vector<Point<dimension - 1>>& z_boundary) = 0;
};
}

namespace Shape {
template <uint dimension>
class Shape {
  protected:
    std::vector<Point<dimension>> nodal_coordinates;

  public:
    Shape(const std::vector<Point<dimension>>& nodal_coordinates) : nodal_coordinates(std::move(nodal_coordinates)) {}

    virtual bool CheckJacobianPositive(const Point<dimension>& point) = 0;

    virtual std::vector<double> GetJdet(const std::vector<Point<dimension>>& points) = 0;
    virtual Array3D<double> GetJinv(const std::vector<Point<dimension>>& points) = 0;
    virtual std::vector<double> GetSurfaceJ(const uint bound_id, const std::vector<Point<dimension>>& points) = 0;
    virtual Array2D<double> GetSurfaceNormal(const uint bound_id, const std::vector<Point<dimension>>& points) = 0;

    virtual std::vector<double> InterpolateNodalValues(const std::vector<double>& nodal_values,
                                                       const std::vector<Point<dimension>>& points) = 0;
    virtual std::vector<Point<dimension>> LocalToGlobalCoordinates(const std::vector<Point<dimension>>& points) = 0;

    virtual void GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) = 0;
};
}

//#define RESL2

#define PI 3.14159265359

#define N_DIV 1                // postproc elem div
#define DEFAULT_ID 4294967295  // max uint as default id
#define INTERNAL 255           // max uchar as default bound type: internal
#define DISTRIBUTED 254

enum GlobalCoord : uchar {
    x = 0,
    y = 1,
    z = 2
};

enum LocalCoordTri : uchar {
    z1 = 0,
    z2 = 1,
    z3 = 2
};

enum LocalCoordQuad : uchar {
    n1 = 0,
    n2 = 1,
    n3 = 2
};

enum VTKElementTypes : uchar {
    straight_triangle = 5
};

#endif
