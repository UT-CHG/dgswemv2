#ifndef IHDG_SWE_BC_OUTFLOW_HPP
#define IHDG_SWE_BC_OUTFLOW_HPP

#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"
#include "problem/SWE/discretization_EHDG/boundary_conditions/compute_bc_trace.hpp"

namespace SWE {
namespace IHDG {
namespace BC {
class Outflow {
  private:
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
    void ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename BoundaryType>
void Outflow::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);

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
void Outflow::ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& bound = edge_bound.boundary;

    auto& state    = bound.data.state[0];
    auto& boundary = bound.data.boundary[bound.bound_id];

    auto& edge_state = edge_bound.edge_data.edge_state;

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        boundary.q_at_gp[var] = bound.ComputeUgp(state.q[var]);
    }

    set_constant(edge_state.q_hat, 0.0);

    // this is not implemented
    std::cout << "outflow BC is not implemented\n";
    abort();
}

template <typename StepperType, typename EdgeBoundaryType>
void Outflow::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    auto& q_hat_at_gp    = edge_internal.q_hat_at_gp;
    auto& aux_hat_at_gp  = edge_internal.aux_hat_at_gp;
    auto& surface_normal = edge_bound.boundary.surface_normal;

    get_Aplus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aplus);
    get_dAplus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dze);
    get_dAplus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqx);
    get_dAplus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAplus_dqy);

    get_Aminus(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->Aminus);
    get_dAminus_dze(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dze);
    get_dAminus_dqx(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqx);
    get_dAminus_dqy(q_hat_at_gp, aux_hat_at_gp, surface_normal, this->dAminus_dqy);

    StatVector<double, SWE::n_variables> q_minus_q_hat;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        auto q_hat = column(q_hat_at_gp, gp);

        StatMatrix<double, SWE::n_variables, SWE::n_variables> dB_dq_hat = this->Aminus[gp] - this->Aplus[gp];

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            q_minus_q_hat[var] = boundary.q_at_gp[var][gp] - q_hat[var];
        }
        auto& q_inf_minus_q_hat = q_minus_q_hat;

        column(dB_dq_hat, SWE::Variables::ze) +=
            this->dAplus_dze[gp] * q_minus_q_hat - this->dAminus_dze[gp] * q_inf_minus_q_hat;
        column(dB_dq_hat, SWE::Variables::qx) +=
            this->dAplus_dqx[gp] * q_minus_q_hat - this->dAminus_dqx[gp] * q_inf_minus_q_hat;
        column(dB_dq_hat, SWE::Variables::qy) +=
            this->dAplus_dqy[gp] * q_minus_q_hat - this->dAminus_dqy[gp] * q_inf_minus_q_hat;

        column(boundary.delta_global_kernel_at_gp, gp)          = flatten<double>(this->Aplus[gp]);
        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = flatten<double>(dB_dq_hat);
        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            this->Aplus[gp] * q_minus_q_hat - this->Aminus[gp] * q_inf_minus_q_hat;
    }
}
}
}
}

#endif
