#ifndef EHDG_SWE_BC_LAND_HPP
#define EHDG_SWE_BC_LAND_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"
#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"

namespace SWE {
namespace EHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {}

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    auto n_x = row(edge_bound.boundary.surface_normal, GlobalCoord::x);
    auto n_y = row(edge_bound.boundary.surface_normal, GlobalCoord::y);

    auto qn = vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qx), n_x) +
              vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qy), n_y);

    row(edge_internal.q_hat_at_gp, SWE::Variables::ze) = row(boundary.q_at_gp, SWE::Variables::ze);
    row(edge_internal.q_hat_at_gp, SWE::Variables::qx) =
        row(boundary.q_at_gp, SWE::Variables::qx) - vec_cw_mult(qn, n_x);
    row(edge_internal.q_hat_at_gp, SWE::Variables::qy) =
        row(boundary.q_at_gp, SWE::Variables::qy) - vec_cw_mult(qn, n_y);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    /* Compute trace flux */

    SWE::get_Fn(edge_internal.q_hat_at_gp,
                edge_internal.aux_hat_at_gp,
                edge_bound.boundary.surface_normal,
                boundary.F_hat_at_gp);

    /* Add stabilization parameter terms */

    SWE::get_tau_LF(
        edge_internal.q_hat_at_gp, edge_internal.aux_hat_at_gp, edge_bound.boundary.surface_normal, edge_internal.tau);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        column(boundary.F_hat_at_gp, gp) +=
            edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
    }
}
}
}
}

#endif