#ifndef IHDG_SWE_PROC_EDGE_BOUND_HPP
#define IHDG_SWE_PROC_EDGE_BOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeBoundaryType>
void Problem::prepare_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    edge_bound.boundary.boundary_condition.ComputeNumericalFlux(edge_bound);

    add_dF_hat_tau_terms_boundary_LF(edge_bound);

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);
}
}
}

#endif