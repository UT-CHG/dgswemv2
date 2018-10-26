#ifndef EHDG_GN_BC_LAND_HPP
#define EHDG_GN_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/Green-Naghdi/discretization_EHDG/stabilization_parameters/ehdg_gn_stabilization_params.hpp"

namespace GN {
namespace EHDG {
namespace BC {
class Land : public SWE::EHDG::BC::Land {
  public:
    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeGlobalKernelsDC(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double nx, ny;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        nx = edge_bound.boundary.surface_normal(GlobalCoord::x, gp);
        ny = edge_bound.boundary.surface_normal(GlobalCoord::y, gp);

        column(edge_internal.w1_hat_w1_hat_kernel_at_gp, gp) = -IdentityVector<double>(GN::n_dimensions);

        column(boundary.w1_hat_w1_kernel_at_gp, gp) = IdentityVector<double>(GN::n_dimensions);

        boundary.w1_hat_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::x, gp) -= nx * nx;
        boundary.w1_hat_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::x + GlobalCoord::y, gp) -= nx * ny;
        boundary.w1_hat_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::x, gp) -= nx * ny;
        boundary.w1_hat_w1_kernel_at_gp(GN::n_dimensions * GlobalCoord::y + GlobalCoord::y, gp) -= ny * ny;
    }
}
}
}
}

#endif