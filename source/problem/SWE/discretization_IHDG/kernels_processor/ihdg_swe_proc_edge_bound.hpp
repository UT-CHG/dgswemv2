#ifndef IHDG_SWE_PROC_EDGE_BOUND_HPP
#define IHDG_SWE_PROC_EDGE_BOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeBoundaryType>
void Problem::local_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& internal = edge_bound.boundary.data.internal;
    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_bound.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp(SWE::Auxiliaries::h, gp) =
            edge_internal.q_hat_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
    }

    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_bound_LF(edge_bound);

    // Add (dtau/dq_hat * del_q - tau) term to dF_hat_dq_hat
    // and tau term to dF_hat_dq
    add_dF_hat_tau_terms_bound_LF(edge_bound);

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); dof_j++) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_bound.boundary.IntegrationPhiPhi(dof_j, dof_i, boundary.dF_hat_dq_at_gp));
        }

        subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -edge_bound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); dof_j++) {
            submatrix(boundary.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_bound.IntegrationPhiLambda(dof_i, dof_j, boundary.dF_hat_dq_hat_at_gp));
        }
    }
}

template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); dof_j++) {
            submatrix(boundary.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_bound.IntegrationPhiLambda(dof_j, dof_i, boundary.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); dof_j++) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables, SWE::n_variables, SWE::n_variables * SWE::n_variables>(
                    edge_bound.IntegrationLambdaLambda(dof_j, dof_i, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }
}
}
}

#endif