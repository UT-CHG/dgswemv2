#ifndef IHDG_SWE_IS_INTERNAL_HPP
#define IHDG_SWE_IS_INTERNAL_HPP

namespace SWE {
namespace IHDG {
namespace ISP {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename EdgeInterfaceType>
    void ComputeInitTrace(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);
};

template <typename InterfaceType>
void Internal::Initialize(InterfaceType& intface) {} /*nothing to initialize*/

template <typename EdgeInterfaceType>
void Internal::ComputeInitTrace(EdgeInterfaceType& edge_int) {
    auto& intface = edge_int.interface;

    auto& state_in    = intface.data_in.state[0];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[0];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_int.interface.surface_normal_in;

    boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
    boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);

    set_constant(edge_state.q_hat, 0.0);

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        edge_internal.q_hat_at_gp = edge_int.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);

        SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);
        SWE::get_dtau_dze_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dze);
        SWE::get_dtau_dqx_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqx);
        SWE::get_dtau_dqy_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqy);

        StatVector<double, SWE::n_variables> del_q;
        StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

        uint gp_ex;
        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

            del_q =
                column(boundary_in.q_at_gp, gp) + column(boundary_ex.q_at_gp, gp_ex) - 2.0 * column(q_hat_at_gp, gp);

            column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
            column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
            column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

            column(edge_internal.rhs_global_kernel_at_gp, gp) = edge_internal.tau[gp] * del_q;

            column(edge_internal.delta_hat_global_kernel_at_gp, gp) =
                flatten<double>(dtau_delq - 2.0 * edge_internal.tau[gp]);
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

        solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

        edge_state.q_hat +=
            reshape<double, SWE::n_variables, SO::ColumnMajor>(edge_internal.rhs_global, edge_int.edge_data.get_ndof());

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-12) {
            break;
        }
    }
}

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(boundary_in.delta_global_kernel_at_gp, gp)    = column(boundary_in.dF_hat_dq_at_gp, gp);
        column(boundary_ex.delta_global_kernel_at_gp, gp_ex) = column(boundary_ex.dF_hat_dq_at_gp, gp_ex);

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = column(boundary_in.dF_hat_dq_hat_at_gp, gp);
        column(edge_internal.delta_hat_global_kernel_at_gp, gp) += column(boundary_ex.dF_hat_dq_hat_at_gp, gp_ex);

        column(edge_internal.rhs_global_kernel_at_gp, gp) = column(boundary_in.F_hat_at_gp, gp);
        column(edge_internal.rhs_global_kernel_at_gp, gp) += column(boundary_ex.F_hat_at_gp, gp_ex);
    }
}
}
}
}

#endif