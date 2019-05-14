#ifndef COMPUTE_BC_TRACE_HPP
#define COMPUTE_BC_TRACE_HPP

namespace SWE {
template <typename EdgeBoundaryType>
void compute_bc_trace(EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary           = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];
    auto& boundary_condition = edge_bound.boundary.boundary_condition;

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    StatVector<double, SWE::n_variables> q;
    uint iter = 0;
    while (iter != 100) {
        edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        get_Aplus(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.Aplus);
        get_dAplus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.dAplus_dze);
        get_dAplus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.dAplus_dqx);
        get_dAplus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.dAplus_dqy);

        get_Aminus(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.Aminus);
        get_dAminus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.dAminus_dze);
        get_dAminus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.dAminus_dqx);
        get_dAminus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, boundary_condition.dAminus_dqy);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            for ( uint var = 0; var < SWE::n_variables; ++var ) {
                q[var] = boundary.q_at_gp[var][gp];
            }
            auto q_hat = column(edge_internal.q_hat_at_gp, gp);
            auto q_inf = column(boundary_condition.q_ex, gp);

            StatMatrix<double, SWE::n_variables, SWE::n_variables> dB_dq_hat =
                boundary_condition.Aminus[gp] - boundary_condition.Aplus[gp];

            column(dB_dq_hat, SWE::Variables::ze) +=
                boundary_condition.dAplus_dze[gp] * (q - q_hat) - boundary_condition.dAminus_dze[gp] * (q_inf - q_hat);
            column(dB_dq_hat, SWE::Variables::qx) +=
                boundary_condition.dAplus_dqx[gp] * (q - q_hat) - boundary_condition.dAminus_dqx[gp] * (q_inf - q_hat);
            column(dB_dq_hat, SWE::Variables::qy) +=
                boundary_condition.dAplus_dqy[gp] * (q - q_hat) - boundary_condition.dAminus_dqy[gp] * (q_inf - q_hat);

            column(edge_internal.delta_hat_global_kernel_at_gp, gp) = flatten<double>(dB_dq_hat);
            column(edge_internal.rhs_global_kernel_at_gp, gp) =
                boundary_condition.Aplus[gp] * (q - q_hat) - boundary_condition.Aminus[gp] * (q_inf - q_hat);
        }

        for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
            for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
                submatrix(edge_internal.delta_hat_global,
                          SWE::n_variables * dof_i,
                          SWE::n_variables * dof_j,
                          SWE::n_variables,
                          SWE::n_variables) =
                    reshape<double, SWE::n_variables>(
                        edge_bound.IntegrationLambdaLambda(dof_i, dof_j, edge_internal.delta_hat_global_kernel_at_gp));
            }

            subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
                -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
        }

        solve_sle(edge_internal.delta_hat_global, edge_internal.rhs_global);

        edge_state.q_hat += reshape<double, SWE::n_variables, SO::ColumnMajor>(edge_internal.rhs_global,
                                                                               edge_bound.edge_data.get_ndof());

        double delta_hat_norm = norm(edge_internal.rhs_global) / edge_internal.rhs_global.size();

        if (delta_hat_norm < 1.0e-12) {
            break;
        }

        ++iter;
    }
}
}

#endif