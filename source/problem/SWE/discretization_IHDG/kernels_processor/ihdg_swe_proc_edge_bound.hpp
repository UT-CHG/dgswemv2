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
        edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_bound_LF(edge_bound);

    // Add (dtau/dq_hat * del_q - tau) term to dF_hat_dq_hat
    // and tau term to dF_hat_dq
    add_dF_hat_tau_terms_bound_LF(edge_bound);

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); dof_j++) {
            blaze::submatrix(internal.delta_local,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) +=
                edge_bound.boundary.IntegrationPhiPhi(dof_i, dof_j, boundary.dF_hat_dq_at_gp);
        }

        blaze::subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -edge_bound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp);
    }

    // Initialize delta_hat container
    boundary.delta_hat_local.resize(SWE::n_variables * edge_bound.boundary.data.get_ndof(),
                                    SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                    false);

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); dof_j++) {
            blaze::submatrix(boundary.delta_hat_local,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) =
                edge_bound.IntegrationPhiLambda(dof_i, dof_j, boundary.dF_hat_dq_hat_at_gp);
        }
    }
}

template <typename EdgeBoundaryType>
void Problem::global_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

    // Initialize delta_hat and rhs containers
    edge_internal.delta_hat_global.resize(
        SWE::n_variables * edge_bound.edge_data.get_ndof(), SWE::n_variables * edge_bound.edge_data.get_ndof(), false);
    edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(), false);

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); dof_j++) {
            blaze::submatrix(edge_internal.delta_hat_global,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) =
                edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp);
        }

        blaze::subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }

    // Initialize delta_hat and rhs containers
    boundary.delta_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                 SWE::n_variables * edge_bound.boundary.data.get_ndof(),
                                 false);

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); dof_j++) {
            blaze::submatrix(boundary.delta_global,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) =
                edge_bound.IntegrationPhiLambda(dof_j, dof_i, boundary.delta_global_kernel_at_gp);
        }
    }
}
}
}

#endif