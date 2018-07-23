#ifndef RKDG_SWE_BC_FLOW_HPP
#define RKDG_SWE_BC_FLOW_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace BC {
class Flow {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;
    DynRowVector<double> qn;

  private:
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<DynRowVector<double>> amplitude;
    std::vector<DynRowVector<double>> phase;

    std::vector<DynRowVector<double>> amplitude_gp;
    std::vector<DynRowVector<double>> phase_gp;

  public:
    Flow() = default;
    Flow(const std::vector<FlowInput>& flow_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void ComputeFlux(const RKStepper& stepper,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     const HybMatrix<double, SWE::n_variables>& q_in,
                     const HybMatrix<double, SWE::n_variables>& aux_in,
                     HybMatrix<double, SWE::n_variables>& F_hat);
};

Flow::Flow(const std::vector<FlowInput>& flow_input) {
    this->frequency    = flow_input[0].frequency;
    this->forcing_fact = flow_input[0].forcing_fact;
    this->equilib_arg  = flow_input[0].equilib_arg;

    uint n_contituents = this->frequency.size();
    uint n_nodes       = flow_input.size();

    this->amplitude.resize(n_contituents);
    this->phase.resize(n_contituents);

    for (uint con = 0; con < n_contituents; con++) {
        this->amplitude[con].resize(n_nodes);
        this->phase[con].resize(n_nodes);

        for (uint node = 0; node < n_nodes; node++) {
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

    for (uint con = 0; con < n_contituents; con++) {
        this->amplitude_gp[con] = bound.ComputeBoundaryNodalUgp(this->amplitude[con]);
        this->phase_gp[con]     = bound.ComputeBoundaryNodalUgp(this->phase[con]);
    }
}

void Flow::ComputeFlux(const RKStepper& stepper,
                       const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                       const HybMatrix<double, SWE::n_variables>& q_in,
                       const HybMatrix<double, SWE::n_variables>& aux_in,
                       HybMatrix<double, SWE::n_variables>& F_hat) {
    // *** //
    std::fill(this->qn.begin(), this->qn.end(), 0.0);

    for (uint con = 0; con < this->frequency.size(); con++) {
        for (uint gp = 0; gp < columns(q_in); ++gp) {
            this->qn[gp] += stepper.GetRamp() * this->forcing_fact[con] * this->amplitude_gp[con][gp] *
                            cos(this->frequency[con] * stepper.GetTimeAtCurrentStage() + this->equilib_arg[con] -
                                this->phase_gp[con][gp]);
        }
    }

    auto n_x = row(surface_normal, GlobalCoord::x);
    auto n_y = row(surface_normal, GlobalCoord::y);

    row(this->q_ex, SWE::Variables::ze) = row(q_in, SWE::Variables::ze);
    row(this->q_ex, SWE::Variables::qx) = cwise_multiplication(qn, n_x);
    row(this->q_ex, SWE::Variables::qy) = cwise_multiplication(qn, n_y);

    for (uint gp = 0; gp < columns(q_in); ++gp) {
        column(F_hat, gp) = LLF_flux(
            Global::g, column(q_in, gp), column(this->q_ex, gp), column(aux_in, gp), column(surface_normal, gp));
    }
}
}
}
}

#endif