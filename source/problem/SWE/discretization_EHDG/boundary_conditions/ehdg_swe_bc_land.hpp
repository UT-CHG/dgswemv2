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

    template <typename StepperType, typename BoundaryType>
    void ComputeBedFlux(const StepperType& stepper, BoundaryType& bound);

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
    if (edge_bound.boundary.data.wet_dry_state.wet) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;
        auto& boundary      = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        const auto n_x = row(edge_bound.boundary.surface_normal, GlobalCoord::x);
        const auto n_y = row(edge_bound.boundary.surface_normal, GlobalCoord::y);
        const auto qn  = vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qx), n_x) +
                        vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::qy), n_y);

        // Assume bath_hat = 0.5*(bath_in +bath_ex)
        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath) = row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) = row(boundary.q_at_gp, SWE::Variables::ze);
        row(edge_internal.q_hat_at_gp, SWE::Variables::qx) =
            row(boundary.q_at_gp, SWE::Variables::qx) - vec_cw_mult(qn, n_x);
        row(edge_internal.q_hat_at_gp, SWE::Variables::qy) =
            row(boundary.q_at_gp, SWE::Variables::qy) - vec_cw_mult(qn, n_y);
        row(edge_internal.q_hat_at_gp, SWE::Variables::hc) = row(boundary.q_at_gp, SWE::Variables::hc);

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
            row(edge_internal.q_hat_at_gp, SWE::Variables::ze) +
            row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

        /* Compute trace flux */
        SWE::get_Fn(edge_internal.q_hat_at_gp,
                    edge_internal.aux_hat_at_gp,
                    edge_bound.boundary.surface_normal,
                    boundary.F_hat_at_gp);

        /* Add stabilization parameter terms */
        SWE::get_tau_LF(edge_internal.q_hat_at_gp,
                        edge_internal.aux_hat_at_gp,
                        edge_bound.boundary.surface_normal,
                        edge_internal.tau);

        for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
            column(boundary.F_hat_at_gp, gp) +=
                edge_internal.tau[gp] * (column(boundary.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
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

    SWE::get_Fn(q_hat, aux_hat, surface_normal, std::move(F_hat));

    SWE::get_tau_LF(q_hat, aux_hat, surface_normal, tau);

    F_hat += tau * (q - q_hat);
}

template <typename StepperType, typename BoundaryType>
void Land::ComputeBedFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];
    set_constant(boundary.qb_hat_at_gp, 0.0);
}
}
}
}

#endif