#ifndef EHDG_SWE_BC_LAND_HPP
#define EHDG_SWE_BC_LAND_HPP

#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"

namespace SWE {
namespace EHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename EdgeBoundaryType>
    void ComputeNumericalFlux(EdgeBoundaryType& edge_bound);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {}

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

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
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

    row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h) =
        row(edge_internal.q_hat_at_gp, SWE::Variables::ze) + row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath);

    /* Compute trace flux */

    auto nx = row(edge_bound.boundary.surface_normal, GlobalCoord::x);
    auto ny = row(edge_bound.boundary.surface_normal, GlobalCoord::y);

    auto u = vec_cw_div(row(edge_internal.q_hat_at_gp, SWE::Variables::qx),
                        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));
    auto v = vec_cw_div(row(edge_internal.q_hat_at_gp, SWE::Variables::qy),
                        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::h));

    auto uuh = vec_cw_mult(u, row(edge_internal.q_hat_at_gp, SWE::Variables::qx));
    auto vvh = vec_cw_mult(v, row(edge_internal.q_hat_at_gp, SWE::Variables::qy));
    auto uvh = vec_cw_mult(u, row(edge_internal.q_hat_at_gp, SWE::Variables::qy));
    auto pe  = Global::g * (0.5 * vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::ze),
                                             row(edge_internal.q_hat_at_gp, SWE::Variables::ze)) +
                           vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::ze),
                                       row(boundary.aux_at_gp, SWE::Auxiliaries::bath)));

    row(boundary.F_hat_at_gp, SWE::Variables::ze) =
        vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::qx), nx) +
        vec_cw_mult(row(edge_internal.q_hat_at_gp, SWE::Variables::qy), ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qx) = vec_cw_mult(uuh + pe, nx) + vec_cw_mult(uvh, ny);
    row(boundary.F_hat_at_gp, SWE::Variables::qy) = vec_cw_mult(uvh, nx) + vec_cw_mult(vvh + pe, ny);

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