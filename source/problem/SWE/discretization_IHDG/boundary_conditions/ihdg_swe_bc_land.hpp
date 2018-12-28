#ifndef IHDG_SWE_BC_LAND_HPP
#define IHDG_SWE_BC_LAND_HPP

namespace SWE {
namespace IHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {} /*nothing to initialize*/

    template <typename EdgeBoundaryType>
    void ComputeInitTrace(EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename EdgeBoundaryType>
void Land::ComputeInitTrace(EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    uint iter = 0;
    while (iter != 100) {
        ++iter;

        q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

        row(aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(q_hat_at_gp, SWE::Variables::ze) + row(aux_hat_at_gp, SWE::Auxiliaries::bath);

        double qn;
        double nx, ny;

        StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            nx = edge_bound.boundary.surface_normal(GlobalCoord::x, gp);
            ny = edge_bound.boundary.surface_normal(GlobalCoord::y, gp);

            qn = boundary.q_at_gp(SWE::Variables::qx, gp) * nx + boundary.q_at_gp(SWE::Variables::qy, gp) * ny;

            column(edge_internal.delta_hat_global_kernel_at_gp, gp) = I_vector;

            column(edge_internal.rhs_global_kernel_at_gp, gp) =
                column(edge_internal.q_hat_at_gp, gp) - column(boundary.q_at_gp, gp);
            edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qx, gp) += qn * nx;
            edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qy, gp) += qn * ny;
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
    }
}

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double qn;
    double nx, ny;

    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        nx = edge_bound.boundary.surface_normal(GlobalCoord::x, gp);
        ny = edge_bound.boundary.surface_normal(GlobalCoord::y, gp);

        qn = boundary.q_at_gp(SWE::Variables::qx, gp) * nx + boundary.q_at_gp(SWE::Variables::qy, gp) * ny;

        column(boundary.delta_global_kernel_at_gp, gp) = -I_vector;

        boundary.delta_global_kernel_at_gp(JacobianVariables::qx_qx, gp) += nx * nx;
        boundary.delta_global_kernel_at_gp(JacobianVariables::qx_qy, gp) += nx * ny;
        boundary.delta_global_kernel_at_gp(JacobianVariables::qy_qx, gp) += nx * ny;
        boundary.delta_global_kernel_at_gp(JacobianVariables::qy_qy, gp) += ny * ny;

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = I_vector;

        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            column(edge_internal.q_hat_at_gp, gp) - column(boundary.q_at_gp, gp);
        edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qx, gp) += qn * nx;
        edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qy, gp) += qn * ny;
    }
}
}
}
}

#endif