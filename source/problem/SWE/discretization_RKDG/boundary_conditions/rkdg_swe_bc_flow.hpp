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
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    Array2D<double> amplitude;
    Array2D<double> phase;

    Array2D<double> amplitude_gp;
    Array2D<double> phase_gp;

  public:
    Flow() = default;
    Flow(const std::vector<FlowInput>& flow_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void ComputeFlux(const RKStepper& stepper,
                     const DynVector<StatVector<double, SWE::n_dimensions>>& surface_normal,
                     const std::vector<StatVector<double, SWE::n_variables>>& q_in,
                     const std::vector<StatVector<double, SWE::n_auxiliaries>>& aux_in,
                     std::vector<StatVector<double, SWE::n_variables>>& F_hat);

    void GetEX(const RKStepper& stepper,
               const uint gp,
               const StatVector<double, SWE::n_dimensions>& surface_normal,
               const StatVector<double, SWE::n_variables>& q_in,
               StatVector<double, SWE::n_variables>& q_ex);
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
    uint n_contituents = this->frequency.size();

    this->amplitude_gp.resize(n_contituents);
    this->phase_gp.resize(n_contituents);

    for (uint con = 0; con < n_contituents; con++) {
        this->amplitude_gp[con].resize(bound.data.get_ngp_boundary(bound.bound_id));
        this->phase_gp[con].resize(bound.data.get_ngp_boundary(bound.bound_id));

        bound.ComputeBoundaryNodalUgp(this->amplitude[con], this->amplitude_gp[con]);
        bound.ComputeBoundaryNodalUgp(this->phase[con], this->phase_gp[con]);
    }
}

void Flow::ComputeFlux(const RKStepper& stepper,
                       const DynVector<StatVector<double, SWE::n_dimensions>>& surface_normal,
                       const std::vector<StatVector<double, SWE::n_variables>>& q_in,
                       const std::vector<StatVector<double, SWE::n_auxiliaries>>& aux_in,
                       std::vector<StatVector<double, SWE::n_variables>>& F_hat) {
    StatVector<double, SWE::n_variables> q_ex;
    for (uint gp = 0; gp < q_in.size(); ++gp) {
        this->GetEX(stepper, gp, surface_normal[gp], q_in[gp], q_ex);

        LLF_flux(Global::g, q_in[gp], q_ex, aux_in[gp], surface_normal[gp], F_hat[gp]);
    }
}

void Flow::GetEX(const RKStepper& stepper,
                 const uint gp,
                 const StatVector<double, SWE::n_dimensions>& surface_normal,
                 const StatVector<double, SWE::n_variables>& q_in,
                 StatVector<double, SWE::n_variables>& q_ex) {
    double qn = 0;
    double qt = 0;

    double frequency;
    double forcing_fact;
    double eq_argument;

    for (uint con = 0; con < this->frequency.size(); con++) {
        frequency    = this->frequency[con];
        forcing_fact = this->forcing_fact[con];
        eq_argument  = this->equilib_arg[con];

        qn += stepper.GetRamp() * forcing_fact * this->amplitude_gp[con][gp] *
              cos(frequency * stepper.GetTimeAtCurrentStage() + eq_argument - this->phase_gp[con][gp]);
    }

    double n_x, n_y, t_x, t_y;

    n_x = surface_normal[GlobalCoord::x];
    n_y = surface_normal[GlobalCoord::y];
    t_x = -n_y;
    t_y = n_x;

    q_ex[SWE::Variables::ze] = q_in[SWE::Variables::ze];
    q_ex[SWE::Variables::qx] = qn * n_x + qt * t_x;
    q_ex[SWE::Variables::qy] = qn * n_y + qt * t_y;
}
}
}
}

#endif