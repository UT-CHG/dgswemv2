#ifndef IHDG_SWE_PROC_EDGE_BOUND_HPP
#define IHDG_SWE_PROC_EDGE_BOUND_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename EdgeBoundaryType>
void Problem::init_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    if (stepper.GetOrder() == 2) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& internal = edge_bound.boundary.data.internal;
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
        auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
        auto& surface_normal = edge_bound.boundary.surface_normal;

        edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        /* Compute fluxes at boundary states */

        auto nx = row(surface_normal, GlobalCoord::x);
        auto ny = row(surface_normal, GlobalCoord::y);

        auto u = vec_cw_div(row(q_hat_at_gp, SWE::Variables::qx), row(aux_hat_at_gp, SWE::Auxiliaries::h));
        auto v = vec_cw_div(row(q_hat_at_gp, SWE::Variables::qy), row(aux_hat_at_gp, SWE::Auxiliaries::h));

        auto uuh = vec_cw_mult(u, row(q_hat_at_gp, SWE::Variables::qx));
        auto vvh = vec_cw_mult(v, row(q_hat_at_gp, SWE::Variables::qy));
        auto uvh = vec_cw_mult(u, row(q_hat_at_gp, SWE::Variables::qy));
        auto pe  = Global::g *
                  (0.5 * vec_cw_mult(row(q_hat_at_gp, SWE::Variables::ze), row(q_hat_at_gp, SWE::Variables::ze)) +
                   vec_cw_mult(row(q_hat_at_gp, SWE::Variables::ze), row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

        // Fn terms
        row(boundary.F_hat_at_gp, SWE::Variables::ze) = vec_cw_mult(row(q_hat_at_gp, SWE::Variables::qx), nx) +
                                                        vec_cw_mult(row(q_hat_at_gp, SWE::Variables::qy), ny);
        row(boundary.F_hat_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
        row(boundary.F_hat_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

        /* Add stabilization parameter terms */

        SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            column(boundary.F_hat_at_gp, gp) +=
                edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
        }

        for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); ++dof_i) {
            subvector(internal.rhs_prev, SWE::n_variables * dof_i, SWE::n_variables) +=
                -edge_bound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp);
        }
    }
}

template <typename StepperType, typename EdgeBoundaryType>
void Problem::local_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& internal = edge_bound.boundary.data.internal;
    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

    /* Compute fluxes at boundary states */

    auto nx = row(surface_normal, GlobalCoord::x);
    auto ny = row(surface_normal, GlobalCoord::y);

    auto u = vec_cw_div(row(q_hat_at_gp, SWE::Variables::qx), row(aux_hat_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(q_hat_at_gp, SWE::Variables::qy), row(aux_hat_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(q_hat_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(q_hat_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(q_hat_at_gp, SWE::Variables::qy));
    auto pe  = Global::g *
              (0.5 * vec_cw_mult(row(q_hat_at_gp, SWE::Variables::ze), row(q_hat_at_gp, SWE::Variables::ze)) +
               vec_cw_mult(row(q_hat_at_gp, SWE::Variables::ze), row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    // Fn terms
    row(boundary.F_hat_at_gp, SWE::Variables::ze) =
        vec_cw_mult(row(q_hat_at_gp, SWE::Variables::qx), nx) + vec_cw_mult(row(q_hat_at_gp, SWE::Variables::qy), ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

    // dFn/dq_hat terms
    set_constant(row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::ze_ze), 0.0);
    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::ze_qx) = nx;
    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::ze_qy) = ny;

    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::qx_ze) =
        vec_cw_mult(-vec_cw_mult(u, u) + Global::g * row(aux_hat_at_gp, SWE::Auxiliaries::h), nx) -
        vec_cw_mult(vec_cw_mult(u, v), ny);
    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::qx_qx) = 2.0 * vec_cw_mult(u, nx) + vec_cw_mult(v, ny);
    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::qx_qy) = vec_cw_mult(u, ny);

    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::qy_ze) =
        -vec_cw_mult(vec_cw_mult(u, v), nx) +
        vec_cw_mult(-vec_cw_mult(v, v) + Global::g * row(aux_hat_at_gp, SWE::Auxiliaries::h), ny);
    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::qy_qx) = vec_cw_mult(v, nx);
    row(boundary.dF_hat_dq_hat_at_gp, JacobianVariables::qy_qy) = vec_cw_mult(u, nx) + 2.0 * vec_cw_mult(v, ny);

    /* Add stabilization parameter terms */

    SWE::get_tau_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.tau);
    SWE::get_dtau_dze_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dze);
    SWE::get_dtau_dqx_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqx);
    SWE::get_dtau_dqy_LF(q_hat_at_gp, aux_hat_at_gp, surface_normal, edge_internal.dtau_dqy);

    StatVector<double, SWE::n_variables> del_q;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        del_q = column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp);

        column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
        column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
        column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

        column(boundary.F_hat_at_gp, gp) += edge_internal.tau[gp] * del_q;
        column(boundary.dF_hat_dq_at_gp, gp) = flatten<double>(edge_internal.tau[gp]);
        column(boundary.dF_hat_dq_hat_at_gp, gp) += flatten<double>(dtau_delq - edge_internal.tau[gp]);
    }

    double theta = stepper.GetTheta();

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(internal.delta_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) +=
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_bound.boundary.IntegrationPhiPhi(dof_j, dof_i, boundary.dF_hat_dq_at_gp));
        }

        subvector(internal.rhs_local, SWE::n_variables * dof_i, SWE::n_variables) +=
            -(1.0 - theta) * edge_bound.boundary.IntegrationPhi(dof_i, boundary.F_hat_at_gp);
    }

    for (uint dof_i = 0; dof_i < edge_bound.boundary.data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
            submatrix(boundary.delta_hat_local,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    (1.0 - theta) * edge_bound.IntegrationPhiLambda(dof_i, dof_j, boundary.dF_hat_dq_hat_at_gp));
        }
    }
}

template <typename StepperType, typename EdgeBoundaryType>
void Problem::global_edge_boundary_kernel(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_bound.boundary.boundary_condition.ComputeGlobalKernels(stepper, edge_bound);

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.boundary.data.get_ndof(); ++dof_j) {
            submatrix(boundary.delta_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_bound.IntegrationPhiLambda(dof_j, dof_i, boundary.delta_global_kernel_at_gp));
        }
    }

    for (uint dof_i = 0; dof_i < edge_bound.edge_data.get_ndof(); ++dof_i) {
        for (uint dof_j = 0; dof_j < edge_bound.edge_data.get_ndof(); ++dof_j) {
            submatrix(edge_internal.delta_hat_global,
                      SWE::n_variables * dof_i,
                      SWE::n_variables * dof_j,
                      SWE::n_variables,
                      SWE::n_variables) =
                reshape<double, SWE::n_variables>(
                    edge_bound.IntegrationLambdaLambda(dof_j, dof_i, edge_internal.delta_hat_global_kernel_at_gp));
        }

        subvector(edge_internal.rhs_global, SWE::n_variables * dof_i, SWE::n_variables) =
            -edge_bound.IntegrationLambda(dof_i, edge_internal.rhs_global_kernel_at_gp);
    }
}
}
}

#endif