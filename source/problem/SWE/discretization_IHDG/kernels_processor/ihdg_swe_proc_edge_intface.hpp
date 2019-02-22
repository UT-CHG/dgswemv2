#ifndef IHDG_SWE_PROC_EDGE_INTFACE_HPP
#define IHDG_SWE_PROC_EDGE_INTFACE_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"

namespace SWE {
namespace IHDG {
template <typename StepperType, typename EdgeInterfaceType>
void Problem::init_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int) {
    if (stepper.GetOrder() == 2) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& internal_in = edge_int.interface.data_in.internal;
        auto& internal_ex = edge_int.interface.data_ex.internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
        auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
        auto& surface_normal = edge_int.interface.surface_normal_in;

        edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        /* Compute fluxes at boundary states */

        SWE::get_Fn(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_in.F_hat_at_gp);

        /* Add stabilization parameter terms */

        SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);

        uint gp_ex;
        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

            column(boundary_ex.F_hat_at_gp, gp_ex) = -column(boundary_in.F_hat_at_gp, gp);

            column(boundary_in.F_hat_at_gp, gp) +=
                edge_internal.tau[gp] * (column(boundary_in.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
            column(boundary_ex.F_hat_at_gp, gp_ex) +=
                edge_internal.tau[gp] * (column(boundary_ex.q_at_gp, gp_ex) - column(edge_internal.q_hat_at_gp, gp));
        }

        for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); ++dof_i) {
            subvector(internal_in.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) +=
                -edge_int.interface.IntegrationPhiIN(dof_i, boundary_in.F_hat_at_gp);
        }

        for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); ++dof_i) {
            subvector(internal_ex.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) +=
                -edge_int.interface.IntegrationPhiEX(dof_i, boundary_ex.F_hat_at_gp);
        }
    }
}

template <typename StepperType, typename EdgeInterfaceType>
void Problem::local_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& internal_in = edge_int.interface.data_in.internal;
    auto& internal_ex = edge_int.interface.data_ex.internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_int.interface.surface_normal_in;

    edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary states */

    SWE::get_Fn(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_in.F_hat_at_gp);

    SWE::get_dFn_dq(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_in.dF_hat_dq_hat_at_gp);

    /* Add stabilization parameter terms */

    SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);
    SWE::get_dtau_dze_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dze);
    SWE::get_dtau_dqx_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqx);
    SWE::get_dtau_dqy_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqy);

    StatVector<double, SWE::n_variables> del_q_in;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq_in;

    StatVector<double, SWE::n_variables> del_q_ex;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq_ex;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        del_q_in = column(boundary_in.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp);
        del_q_ex = column(boundary_ex.q_at_gp, gp_ex) - column(edge_internal.q_hat_at_gp, gp);

        column(dtau_delq_in, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q_in;
        column(dtau_delq_in, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q_in;
        column(dtau_delq_in, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q_in;

        column(dtau_delq_ex, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q_ex;
        column(dtau_delq_ex, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q_ex;
        column(dtau_delq_ex, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q_ex;

        column(boundary_ex.F_hat_at_gp, gp_ex)         = -column(boundary_in.F_hat_at_gp, gp);
        column(boundary_ex.dF_hat_dq_hat_at_gp, gp_ex) = -column(boundary_in.dF_hat_dq_hat_at_gp, gp);

        column(boundary_in.F_hat_at_gp, gp) += edge_internal.tau[gp] * del_q_in;
        column(boundary_in.dF_hat_dq_at_gp, gp) = flatten<double>(edge_internal.tau[gp]);
        column(boundary_in.dF_hat_dq_hat_at_gp, gp) += flatten<double>(dtau_delq_in - edge_internal.tau[gp]);

        column(boundary_ex.F_hat_at_gp, gp_ex) += edge_internal.tau[gp] * del_q_ex;
        column(boundary_ex.dF_hat_dq_at_gp, gp_ex) = flatten<double>(edge_internal.tau[gp]);
        column(boundary_ex.dF_hat_dq_hat_at_gp, gp_ex) += flatten<double>(dtau_delq_ex - edge_internal.tau[gp]);
    }

    double theta = stepper.GetTheta();

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); ++dof_j) {
            submatrix(internal_in.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_int.interface.IntegrationPhiPhiIN(dof_i, dof_j, boundary_in.dF_hat_dq_at_gp));
        }

        subvector(internal_in.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -(1.0 - theta) * edge_int.interface.IntegrationPhiIN(dof_i, boundary_in.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); ++dof_j) {
            submatrix(internal_ex.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_int.interface.IntegrationPhiPhiEX(dof_i, dof_j, boundary_ex.dF_hat_dq_at_gp));
        }

        subvector(internal_ex.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -(1.0 - theta) * edge_int.interface.IntegrationPhiEX(dof_i, boundary_ex.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_in.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            submatrix(boundary_in.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_int.IntegrationPhiLambdaIN(dof_i, dof_j, boundary_in.dF_hat_dq_hat_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.interface.data_ex.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            submatrix(boundary_ex.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_int.IntegrationPhiLambdaEX(dof_i, dof_j, boundary_ex.dF_hat_dq_hat_at_gp));
        }
    }
}

template <typename StepperType, typename EdgeInterfaceType>
void Problem::global_edge_interface_kernel(const StepperType& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    edge_int.interface.specialization.ComputeGlobalKernels(edge_int);

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_in.get_ndof(); ++dof_j) {
            submatrix(boundary_in.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_int.IntegrationPhiLambdaIN(dof_j, dof_i, boundary_in.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.interface.data_ex.get_ndof(); ++dof_j) {
            submatrix(boundary_ex.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_int.IntegrationPhiLambdaEX(dof_j, dof_i, boundary_ex.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_int.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_int.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_int.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_int.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }
}
}
}

#endif