#ifndef EHDG_SWE_BC_FUNCTION_HPP
#define EHDG_SWE_BC_FUNCTION_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"
#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"
#include "compute_bc_trace.hpp"

namespace SWE {
namespace EHDG {
namespace BC {
class Function {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;

    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> Aplus;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAplus_dqy;

    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> Aminus;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dAminus_dqy;

    template <typename EdgeBoundaryType>
    friend void SWE::compute_bc_trace(EdgeBoundaryType& edge_bound);

  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename BoundaryType>
void Function::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);

    this->q_ex.resize(SWE::n_variables, ngp);

    this->Aplus.resize(ngp);
    this->dAplus_dze.resize(ngp);
    this->dAplus_dqx.resize(ngp);
    this->dAplus_dqy.resize(ngp);

    this->Aminus.resize(ngp);
    this->dAminus_dze.resize(ngp);
    this->dAminus_dqx.resize(ngp);
    this->dAminus_dqy.resize(ngp);
}

template <typename StepperType, typename EdgeBoundaryType>
void Function::ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double t = stepper.GetTimeAtCurrentStage();

    this->q_ex = edge_bound.boundary.ComputeFgp([t](Point<2>& pt) {
        double ze = 0.0;
        double qx = 0.0;
        double qy = 0.0;

        Utilities::ignore(ze, qx, qy);

        if (t <= 3.0) {
            ze = cos(PI * t) - 1.0;
        } else {
            ze = -2.0;
        }

        // StatVector<double, SWE::n_variables> q{ze, qx, qy};
        StatVector<double, SWE::n_variables> q(SWE::ic_q(t, pt));

        return q;
    });

    SWE::compute_bc_trace(edge_bound);

    edge_internal.q_hat_at_gp = edge_bound.ComputeUgp(edge_state.q_hat);

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
