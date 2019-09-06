#ifndef EHDG_SWE_BC_FLOW_HPP
#define EHDG_SWE_BC_FLOW_HPP

#include "problem/SWE/problem_flux/swe_flux.hpp"
#include "problem/SWE/problem_jacobian/swe_jacobian.hpp"
#include "compute_bc_trace.hpp"

namespace SWE {
namespace EHDG {
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
    void ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound);
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
void Flow::ComputeNumericalFlux(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state    = edge_bound.edge_data.edge_state;
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    set_constant(this->qn, 0.0);

    for (uint con = 0; con < this->frequency.size(); ++con) {
        for (uint gp = 0; gp < columns(boundary.q_at_gp); ++gp) {
            this->qn[gp] += stepper.GetRamp() * this->forcing_fact[con] * this->amplitude_gp[con][gp] *
                            cos(this->frequency[con] * stepper.GetTimeAtCurrentStage() +
                                (this->equilib_arg[con] - this->phase_gp[con][gp]) * PI / 180);
        }
    }

    auto n_x = row(edge_bound.boundary.surface_normal, GlobalCoord::x);
    auto n_y = row(edge_bound.boundary.surface_normal, GlobalCoord::y);

    row(this->q_ex, SWE::Variables::ze) = row(boundary.q_at_gp, SWE::Variables::ze);
    row(this->q_ex, SWE::Variables::qx) = vec_cw_mult(qn, n_x);
    row(this->q_ex, SWE::Variables::qy) = vec_cw_mult(qn, n_y);

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