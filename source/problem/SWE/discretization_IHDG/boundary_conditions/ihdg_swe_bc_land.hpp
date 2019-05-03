#ifndef IHDG_SWE_BC_LAND_HPP
#define IHDG_SWE_BC_LAND_HPP

#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"

namespace SWE {
namespace IHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {}

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& bound = edge_bound.boundary;

    auto& state    = bound.data.state[0];
    auto& boundary = bound.data.boundary[bound.bound_id];

    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    boundary.q_at_gp = bound.ComputeUgp(state.q);

    auto n_x = row(edge_bound.boundary.surface_normal, GlobalCoord::x);
    auto n_y = row(edge_bound.boundary.surface_normal, GlobalCoord::y);

    auto qn = vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qx), n_x) +
              vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qy), n_y);

    row(edge_internal.q_hat_at_gp, SWE::Variables::ze) = row(boundary.q_at_gp, SWE::Variables::ze);
    row(edge_internal.q_hat_at_gp, SWE::Variables::qx) =
        row(boundary.q_at_gp, SWE::Variables::qx) - vec_cw_mult(qn, n_x);
    row(edge_internal.q_hat_at_gp, SWE::Variables::qy) =
        row(boundary.q_at_gp, SWE::Variables::qy) - vec_cw_mult(qn, n_y);

    edge_state.q_hat = edge_bound.L2Projection(edge_internal.q_hat_at_gp);
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