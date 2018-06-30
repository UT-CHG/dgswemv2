#ifndef EHDG_SWE_BC_LAND_HPP
#define EHDG_SWE_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_EHDG/numerical_fluxes/ehdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace EHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {} /*nothing to initialize*/

    template <typename EdgeBoundaryType>
    void AddDeltaKernelBTerms(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeBoundaryType>
    void AddRHSKernelBTerms(const RKStepper& stepper, EdgeBoundaryType& edge_bound);
};

template <typename EdgeBoundaryType>
void Land::AddDeltaKernelBTerms(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_global = edge_bound.edge_data.edge_global;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_global.delta_ze_hat_kernel_at_gp[Variables::ze][gp] = -1.0;
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qx][gp] = 0.0;
        edge_global.delta_ze_hat_kernel_at_gp[Variables::qy][gp] = 0.0;

        edge_global.delta_qx_hat_kernel_at_gp[Variables::ze][gp] = 0.0;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qx][gp] = -1.0;
        edge_global.delta_qx_hat_kernel_at_gp[Variables::qy][gp] = 0.0;

        edge_global.delta_qy_hat_kernel_at_gp[Variables::ze][gp] = 0.0;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qx][gp] = 0.0;
        edge_global.delta_qy_hat_kernel_at_gp[Variables::qy][gp] = -1.0;
    }
}

template <typename EdgeBoundaryType>
void Land::AddRHSKernelBTerms(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state  = edge_bound.edge_data.edge_state;
    auto& edge_global = edge_bound.edge_data.edge_global;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double qn;
    double nx, ny;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        qn = boundary.qx_at_gp[gp] * nx + boundary.qy_at_gp[gp] * ny;

        edge_global.ze_rhs_kernel_at_gp[gp] = edge_state.ze_hat_at_gp[gp] - boundary.ze_at_gp[gp];
        edge_global.qx_rhs_kernel_at_gp[gp] = edge_state.qx_hat_at_gp[gp] - boundary.qx_at_gp[gp] + qn * nx;
        edge_global.qy_rhs_kernel_at_gp[gp] = edge_state.qy_hat_at_gp[gp] - boundary.qy_at_gp[gp] + qn * ny;
    }
}
}
}
}

#endif