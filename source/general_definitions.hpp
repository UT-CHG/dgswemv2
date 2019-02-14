#ifndef GENERAL_DEFINITIONS
#define GENERAL_DEFINITIONS

#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <assert.h>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>
#include <stdlib.h>
#include <time.h>

using uint  = unsigned int;
using uchar = unsigned char;

template <uint dimension>
using Point = std::array<double, dimension>;
template <typename type>
using Array2D = std::vector<std::vector<type>>;
template <typename type>
using Array3D = std::vector<std::vector<std::vector<type>>>;
template <typename type>
using Array4D = std::vector<std::vector<std::vector<std::vector<type>>>>;

#include "utilities/linear_algebra.hpp"
#include "utilities/edge_types.hpp"

#ifdef HAS_HPX
//#include "simulation/hpx/load_balancer/serialization_headers.hpp"
#endif

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
    virtual DynMatrix<double> GetPhi(const uint p, const std::vector<Point<dimension>>& points) = 0;
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
    virtual StatVector<DynMatrix<double>, dimension> GetDPhi(const uint p,
							     const std::vector<Point<dimension>>& points) = 0;

    /**
     * Obtain the inverted mass matrix of basis functions of polynomial order p.
     *
     * @param p polynomial order
     * @param return a pair with a boolean, which states whether the basis is orthogonal,
     *        and a 2-dimensional array corresponding to the mass matrix over the master element
     */
    virtual DynMatrix<double> GetMinv(const uint p) = 0;

    virtual DynRowVector<double> ProjectBasisToLinear(const DynRowVector<double>& u) =0;
    virtual DynRowVector<double> ProjectLinearToBasis(const uint ndof, const DynRowVector<double>& u_lin) = 0;

    virtual DynMatrix<double> ProjectBasisToLinear(const DynMatrix<double>& u)                      = 0;
    virtual DynMatrix<double> ProjectLinearToBasis(const uint ndof, const DynMatrix<double>& u_lin) = 0;
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
     * Returns a vector of weights and quadrature points.
     *
     * @param p Polynomial order for which the rule should return exact results.
     * @return Number of Gauss points
     */
    virtual std::pair<DynVector<double>, std::vector<Point<dimension>>> GetRule(const uint p) = 0;

    /**
     * Returns the number of Gauss points required for a rule of strength p
     * GetNumGP returns the number of quadrature points for a rule of strength p without
     * without having to construct the rule.
     *
     * @param p Polynomial order for which the rule should return exact results.
     * @return Pair of weights and quadrature points on the master triangle.
     */
    virtual uint GetNumGP(const uint p) = 0;
};
}

namespace Master {
template <uint dimension>
class Master {
  public:
    uint p;

    uint ndof;
    uint ngp;
    uint nvrtx;
    uint nbound;

    std::pair<DynVector<double>, std::vector<Point<dimension>>> integration_rule;

    DynVector<double> chi_baryctr;
    DynMatrix<double> chi_midpts;

    DynMatrix<double> chi_gp;
    DynMatrix<double, SO::ColumnMajor> phi_gp;

    std::array<DynMatrix<double>, dimension> dchi_gp;
    StatVector<DynMatrix<double>, dimension> dphi_gp;

    DynMatrix<double> int_phi_fact;
    DynMatrix<double> int_phi_phi_fact;
    StatVector<DynMatrix<double>, dimension> int_dphi_fact;

    DynMatrix<double, SO::ColumnMajor> m_inv;

    DynMatrix<double> phi_postprocessor_cell;
    DynMatrix<double> phi_postprocessor_point;

  public:
    Master() = default;
    Master(const uint p) : p(p) {}

    virtual std::vector<Point<dimension>> BoundaryToMasterCoordinates(
        const uint bound_id,
        const std::vector<Point<dimension - 1>>& z_boundary) = 0;

    // The following member methods have to be defined (cannot define templated member methods as virual)
    // template <typename InputArrayType>
    // decltype(auto) ComputeLinearUbaryctr(const InputArrayType& u_lin);
    // template <typename InputArrayType>
    // decltype(auto) ComputeLinearUmidpts(const InputArrayType& u_lin);
    // template <typename InputArrayType>
    // decltype(auto) ComputeLinearUvrtx(const InputArrayType& u_lin);
};
}

namespace Shape {
template <uint dimension>
class Shape {
  public:
    std::vector<Point<3>> nodal_coordinates;

    DynMatrix<double> psi_gp;
    std::array<DynMatrix<double>, dimension> dpsi_gp;

  public:
    Shape() = default;
    Shape(std::vector<Point<3>>&& nodal_coordinates) : nodal_coordinates(std::move(nodal_coordinates)) {}

    virtual ~Shape() = default;

    virtual std::vector<uint> GetBoundaryNodeID(const uint bound_id, const std::vector<uint> node_ID) = 0;

    virtual Point<dimension> GetBarycentricCoordinates()           = 0;
    virtual std::vector<Point<dimension>> GetMidpointCoordinates() = 0;

    virtual DynVector<double> GetJdet(const std::vector<Point<dimension>>& points)                          = 0;
    virtual DynVector<double> GetSurfaceJ(const uint bound_id, const std::vector<Point<dimension>>& points) = 0;

    virtual AlignedVector<StatMatrix<double, dimension, dimension>> GetJinv(
        const std::vector<Point<dimension>>& points) = 0;

    virtual AlignedVector<StatVector<double, dimension>> GetSurfaceNormal(
        const uint bound_id,
        const std::vector<Point<dimension>>& points) = 0;

    virtual DynMatrix<double> GetPsi(const std::vector<Point<dimension>>& points)                         = 0;
    virtual std::array<DynMatrix<double>, dimension> GetDPsi(const std::vector<Point<dimension>>& points) = 0;

    virtual DynMatrix<double> GetBoundaryPsi(const uint bound_id, const std::vector<Point<dimension - 1>>& points) = 0;

    virtual std::vector<Point<dimension>> LocalToGlobalCoordinates(const std::vector<Point<dimension>>& points) = 0;

    virtual void GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) = 0;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & nodal_coordinates;
        // clang-format on
    }
    HPX_SERIALIZATION_POLYMORPHIC_ABSTRACT(Shape);
#endif
};
}

/**
 * Base class for the stepper types for a simulation.
 */
class Stepper {
  public:
    /**
     * Get the order accuracy of the stepper
     */
    virtual uint GetOrder() const = 0;
    /**
     * Get the number of stages per timestep
     */
    virtual uint GetNumStages() const = 0;
    /**
     * Get the size of the timestep (in seconds)
     */
    virtual double GetDT() const = 0;

    /**
     * Set the size of the timestep (in seconds)
     */
    virtual void SetDT(double dt) = 0;

    /**
     * Get the current step number
     */
    virtual uint GetStep() const = 0;
    /**
     * Get the current stage number
     */
    virtual uint GetStage() const = 0;
    /**
     * Get the current timestamp
     * The timestamp is the total number of stages that have been executed up to this point.
     */
    virtual uint GetTimestamp() const = 0;

    /**
     * Get the simulated time at the current step and stage
     */
    virtual double GetTimeAtCurrentStage() const = 0;
    /**
     * Get ramp factor
     * Often for stability reasons, we scale boundary or source terms by a number that goes from 0 to 1
     * as the simulation begins. We refer to this factor as ramp. In this stepper, we
     * use a hyperbolic tangent function. When `t = ramp_duration`, `ramp = tanh(2)`.
     */
    virtual double GetRamp() const = 0;

    /**
     * Prefix incrementor advances the stepper by one stage
     * This operation will update all of the internal states of the stepper by one stage.
     */
    virtual Stepper& operator++() = 0;

    // The following member method has to be defined (cannot define templated member methods as virual)
    /**
     * This operation will do one time stage update for an element
     */
    // template <typename ElementType>
    // void UpdateState(ElementType& elt) const;
};

#define PI 3.14159265359

#define N_DIV 2                // postproc elem div
#define DEFAULT_ID 4294967295  // max uint as default id

enum CoordinateSystem : uchar { cartesian = 0, polar = 1, spherical = 2 };

enum GlobalCoord : uchar { x = 0, y = 1, z = 2 };

enum LocalCoordLin : uchar { l1 = 0, l2 = 1, l3 = 2 };

enum LocalCoordTri : uchar { z1 = 0, z2 = 1, z3 = 2 };

enum LocalCoordQuad : uchar { n1 = 0, n2 = 1, n3 = 2 };

enum VTKElementTypes : uchar { straight_triangle = 5 };

#endif
