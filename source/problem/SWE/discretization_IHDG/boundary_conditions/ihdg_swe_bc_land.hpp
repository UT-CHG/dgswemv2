#ifndef IHDG_SWE_BC_LAND_HPP
#define IHDG_SWE_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"

namespace SWE {
namespace IHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {} /*nothing to initialize*/

    template <typename EdgeBoundaryType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeBoundaryType& edge_bound);
};

template <typename EdgeBoundaryType>
void Land::ComputeGlobalKernels(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double qn;
    double nx, ny;

    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        qn = boundary.q_at_gp(SWE::Variables::qx, gp) * nx + boundary.q_at_gp(SWE::Variables::qy, gp) * ny;

        boundary.delta_global_kernel_at_gp[gp] = -I_vector;

        boundary.delta_global_kernel_at_gp[gp][JacobianVariables::qx_qx] += nx * nx;
        boundary.delta_global_kernel_at_gp[gp][JacobianVariables::qx_qy] += nx * ny;
        boundary.delta_global_kernel_at_gp[gp][JacobianVariables::qy_qx] += nx * ny;
        boundary.delta_global_kernel_at_gp[gp][JacobianVariables::qy_qy] += ny * ny;

        row(edge_internal.delta_hat_global_kernel_at_gp, gp) = I_vector;

        row(edge_internal.rhs_global_kernel_at_gp, gp) = row(edge_internal.q_hat_at_gp, gp) - boundary.q_at_gp[gp];
        edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qx, gp) += qn * nx;
        edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qy, gp) += qn * ny;
    }
}
}
}
}

#endif