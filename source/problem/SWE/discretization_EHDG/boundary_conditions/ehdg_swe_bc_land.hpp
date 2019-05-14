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

    void ComputeFlux(const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                     const Column<HybMatrix<double, SWE::n_variables>>& q,
                     const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_hat,
                     StatMatrix<double, SWE::n_variables, SWE::n_variables>& tau,
                     Column<HybMatrix<double, SWE::n_variables>>&& q_hat,
                     Column<HybMatrix<double, SWE::n_variables>>&& F_hat);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {}

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    auto n_x = edge_bound.boundary.surface_normal[GlobalCoord::x];
    auto n_y = edge_bound.boundary.surface_normal[GlobalCoord::y];

    auto qn = vec_cw_mult(boundary.q_at_gp[SWE::Variables::qx], n_x) +
              vec_cw_mult(boundary.q_at_gp[SWE::Variables::qy], n_y);

    row(edge_internal.q_hat_at_gp, SWE::Variables::ze) = boundary.q_at_gp[SWE::Variables::ze];
    row(edge_internal.q_hat_at_gp, SWE::Variables::qx) =
        boundary.q_at_gp[SWE::Variables::qx] - vec_cw_mult(qn, n_x);
    row(edge_internal.q_hat_at_gp, SWE::Variables::qy) =
        boundary.q_at_gp[SWE::Variables::qy] - vec_cw_mult(qn, n_y);

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

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            for ( uint w = 0; w < SWE::n_variables; ++w ) {
                boundary.F_hat_at_gp[var][gp] +=
                    edge_internal.tau[gp](var, w) * (boundary.q_at_gp[w][gp] - edge_internal.q_hat_at_gp(w, gp));
            }
        }
    }
}

void Land::ComputeFlux(const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                       const Column<HybMatrix<double, SWE::n_variables>>& q,
                       const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_hat,
                       StatMatrix<double, SWE::n_variables, SWE::n_variables>& tau,
                       Column<HybMatrix<double, SWE::n_variables>>&& q_hat,
                       Column<HybMatrix<double, SWE::n_variables>>&& F_hat) {
    // *** //
    double n_x = surface_normal[GlobalCoord::x];
    double n_y = surface_normal[GlobalCoord::y];

    double qn = q[SWE::Variables::qx] * n_x + q[SWE::Variables::qy] * n_y;

    q_hat[SWE::Variables::ze] = q[SWE::Variables::ze];
    q_hat[SWE::Variables::qx] = q[SWE::Variables::qx] - qn * n_x;
    q_hat[SWE::Variables::qy] = q[SWE::Variables::qy] - qn * n_y;

    SWE::get_Fn(Global::g, q_hat, aux_hat, surface_normal, std::move(F_hat));

    SWE::get_tau_LF(Global::g, q_hat, aux_hat, surface_normal, tau);

    F_hat += tau * (q - q_hat);
}
}
}
}

#endif