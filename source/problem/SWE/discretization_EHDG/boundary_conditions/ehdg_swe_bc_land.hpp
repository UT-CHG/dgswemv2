#ifndef EHDG_SWE_BC_LAND_HPP
#define EHDG_SWE_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_EHDG/stabilization_parameters/ehdg_swe_stabilization_params.hpp"

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

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

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
}

template <typename EdgeBoundaryType>
void Land::ComputeNumericalFlux(EdgeBoundaryType& edge_bound) {
    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    boundary.F_hat_at_gp = boundary.Fn_at_gp;

    add_F_hat_tau_terms_bound_LF(edge_bound);
}
}
}
}

#endif