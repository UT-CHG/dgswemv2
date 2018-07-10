#ifndef IHDG_SWE_PROC_EDGE_BOUND_HPP
#define IHDG_SWE_PROC_EDGE_BOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeBoundaryType>
void Problem::prepare_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_bound_LF(edge_bound);

    // Add dtau/dq * del_q - tau * del_q term to dF_hat_dq_hat
    // and tau * del_q term to dF_hat_dq
    add_dF_hat_tau_terms_bound_LF(edge_bound);

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);
}
}
}

#endif