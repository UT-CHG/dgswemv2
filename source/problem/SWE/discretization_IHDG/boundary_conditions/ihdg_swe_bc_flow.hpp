#ifndef IHDG_SWE_BC_FLOW_HPP
#define IHDG_SWE_BC_FLOW_HPP

#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"
#include "problem/SWE/discretization_EHDG/boundary_conditions/compute_bc_trace.hpp"

namespace SWE {
namespace IHDG {
namespace BC {
class Flow {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;
    DynRowVector<double> qn;

    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<DynRowVector<double>> amplitude;
    std::vector<DynRowVector<double>> phase;

    std::vector<DynRowVector<double>> amplitude_gp;
    std::vector<DynRowVector<double>> phase_gp;

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
    Flow() = default;
    Flow(const std::vector<FlowNode>& flow_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound);

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

Flow::Flow(const std::vector<FlowNode>& flow_input) {
    this->frequency    = flow_input[0].frequency;
    this->forcing_fact = flow_input[0].forcing_fact;
    this->equilib_arg  = flow_input[0].equilib_arg;

    uint n_contituents = this->frequency.size();
    uint n_nodes       = flow_input.size();

    this->amplitude.resize(n_contituents);
    this->phase.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude[con].resize(n_nodes);
        this->phase[con].resize(n_nodes);

        for (uint node = 0; node < n_nodes; ++node) {
            this->amplitude[con][node] = flow_input[node].amplitude[con];
            this->phase[con][node]     = flow_input[node].phase[con];
        }
    }
}

template <typename BoundaryType>
void Flow::Initialize(BoundaryType& bound) {
    uint ngp           = bound.data.get_ngp_boundary(bound.bound_id);
    uint n_contituents = this->frequency.size();

    this->q_ex.resize(SWE::n_variables, ngp);
    this->qn.resize(ngp);

    this->amplitude_gp.resize(n_contituents);
    this->phase_gp.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude_gp[con] = bound.ComputeBoundaryNodalUgp(this->amplitude[con]);
        this->phase_gp[con]     = bound.ComputeBoundaryNodalUgp(this->phase[con]);
    }

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
void Flow::ComputeInitTrace(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& bound = edge_bound.boundary;

    auto& state    = bound.data.state[0];
    auto& boundary = bound.data.boundary[bound.bound_id];

    auto& edge_state = edge_bound.edge_data.edge_state;

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        boundary.q_at_gp[var] = bound.ComputeUgp(state.q[var]);
    }

    set_constant(edge_state.q_hat, 0.0);

    set_constant(this->qn, 0.0);

    for (uint con = 0; con < this->frequency.size(); ++con) {
        for (uint gp = 0; gp < boundary.q_at_gp[0].size(); ++gp) {
            this->qn[gp] += stepper.GetRamp() * this->forcing_fact[con] * this->amplitude_gp[con][gp] *
                            cos(this->frequency[con] * stepper.GetTimeAtCurrentStage() + this->equilib_arg[con] -
                                this->phase_gp[con][gp]);
        }
    }

    const auto& n_x = edge_bound.boundary.surface_normal[GlobalCoord::x];
    const auto& n_y = edge_bound.boundary.surface_normal[GlobalCoord::y];

    row(this->q_ex, SWE::Variables::ze) = boundary.q_at_gp[SWE::Variables::ze];
    row(this->q_ex, SWE::Variables::qx) = vec_cw_mult(qn, n_x);
    row(this->q_ex, SWE::Variables::qy) = vec_cw_mult(qn, n_y);

    SWE::compute_bc_trace(edge_bound);
}

template <typename StepperType, typename EdgeBoundaryType>
void Flow::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
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

    set_constant(this->qn, 0.0);

    for (uint con = 0; con < this->frequency.size(); ++con) {
        for (uint gp = 0; gp < boundary.q_at_gp[0].size(); ++gp) {
            this->qn[gp] += stepper.GetRampNext() * this->forcing_fact[con] * this->amplitude_gp[con][gp] *
                            cos(this->frequency[con] * stepper.GetTimeAtNextStage() + this->equilib_arg[con] -
                                this->phase_gp[con][gp]);
        }
    }

    const auto& n_x = surface_normal[GlobalCoord::x];
    const auto& n_y = surface_normal[GlobalCoord::y];

    row(this->q_ex, SWE::Variables::ze) = boundary.q_at_gp[SWE::Variables::ze];
    row(this->q_ex, SWE::Variables::qx) = vec_cw_mult(qn, n_x);
    row(this->q_ex, SWE::Variables::qy) = vec_cw_mult(qn, n_y);

    StatVector<double, SWE::n_variables> q_minus_q_hat;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        auto q_hat = column(q_hat_at_gp, gp);
        auto q_inf = column(this->q_ex, gp);

        StatMatrix<double, SWE::n_variables, SWE::n_variables> dB_dq_hat = this->Aminus[gp] - this->Aplus[gp];

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            q_minus_q_hat[var] = boundary.q_at_gp[var][gp] - q_hat[var];
        }

        column(dB_dq_hat, SWE::Variables::ze) +=
            this->dAplus_dze[gp] * q_minus_q_hat - this->dAminus_dze[gp] * (q_inf - q_hat);
        column(dB_dq_hat, SWE::Variables::qx) +=
            this->dAplus_dqx[gp] * q_minus_q_hat - this->dAminus_dqx[gp] * (q_inf - q_hat);
        column(dB_dq_hat, SWE::Variables::qy) +=
            this->dAplus_dqy[gp] * q_minus_q_hat - this->dAminus_dqy[gp] * (q_inf - q_hat);

        column(boundary.delta_global_kernel_at_gp, gp)          = flatten<double>(this->Aplus[gp]);
        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = flatten<double>(dB_dq_hat);
        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            this->Aplus[gp] * q_minus_q_hat - this->Aminus[gp] * (q_inf - q_hat);
    }
}
}
}
}

#endif