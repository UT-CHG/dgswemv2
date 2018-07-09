#ifndef EHDG_SWE_BC_LAND_HPP
#define EHDG_SWE_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"

namespace SWE {
namespace EHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {} /*nothing to initialize*/

    template <typename EdgeBoundaryType>
    void ComputeGlobalKernels(const RKStepper& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeBoundaryType>
    void ComputeNumericalFlux(EdgeBoundaryType& edge_bound);
};

template <typename EdgeBoundaryType>
void Land::ComputeGlobalKernels(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;
    auto& edge_global   = edge_bound.edge_data.edge_global;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_global.delta_hat_kernel_at_gp[gp] = -blaze::IdentityMatrix<double>(SWE::n_variables);
    }

    double qn;
    double nx, ny;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        qn = boundary.q_at_gp[gp][SWE::Variables::qx] * nx + boundary.q_at_gp[gp][SWE::Variables::qy] * ny;

        edge_global.rhs_kernel_at_gp[gp][SWE::Variables::ze] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::ze] - boundary.q_at_gp[gp][SWE::Variables::ze];
        edge_global.rhs_kernel_at_gp[gp][SWE::Variables::qx] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] - boundary.q_at_gp[gp][SWE::Variables::qx] + qn * nx;
        edge_global.rhs_kernel_at_gp[gp][SWE::Variables::qy] =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] - boundary.q_at_gp[gp][SWE::Variables::qy] + qn * ny;
    }
}

template <typename EdgeBoundaryType>
void Land::ComputeNumericalFlux(EdgeBoundaryType& edge_bound) {
    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    boundary.F_hat_at_gp = boundary.Fn_at_gp;

    // Add tau terms
    add_flux_tau_terms_bound_LF(edge_bound);
}
}
}
}

#endif