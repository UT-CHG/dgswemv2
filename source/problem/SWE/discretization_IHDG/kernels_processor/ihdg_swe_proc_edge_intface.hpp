#ifndef IHDG_SWE_PROC_EDGE_INTFACE_HPP
#define IHDG_SWE_PROC_EDGE_INTFACE_HPP

#include <eigen3/Eigen/Dense>

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
void Problem::local_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& local_in = edge_int.interface.data_in.local;
    auto& local_ex = edge_int.interface.data_ex.local;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.ComputeUgp(edge_state.q_hat, edge_internal.q_hat_at_gp);

    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] + boundary_in.aux_at_gp[gp][SWE::Auxiliaries::bath];
    }

    // Add tau * del_q terms to F_hat
    add_F_hat_tau_terms_intface_LF(edge_int);

    // Add dtau/dq * del_q - tau * del_q term to dF_hat_dq_hat
    // and tau * del_q term to dF_hat_dq
    add_dF_hat_tau_terms_intface_LF(edge_int);

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); dof_j++) {
            blaze::submatrix(local_in.delta_matrix,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) +=
                edge_int.interface.IntegrationPhiPhiIN(dof_j, dof_i, boundary_in.dF_hat_dq_at_gp);
        }

        blaze::subvector(local_in.rhs, SWE::n_variables * dof_i, SWE::n_variables) +=
            -edge_int.interface.IntegrationPhiIN(dof_i, boundary_in.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); dof_i++) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); dof_j++) {
            blaze::submatrix(local_ex.delta_matrix,
                             SWE::n_variables * dof_i,
                             SWE::n_variables * dof_j,
                             SWE::n_variables,
                             SWE::n_variables) +=
                edge_int.interface.IntegrationPhiPhiEX(dof_j, dof_i, boundary_ex.dF_hat_dq_at_gp);
        }

        blaze::subvector(local_ex.rhs, SWE::n_variables * dof_i, SWE::n_variables) +=
            -edge_int.interface.IntegrationPhiEX(dof_i, boundary_ex.F_hat_at_gp);
    }
}
}
}

#endif