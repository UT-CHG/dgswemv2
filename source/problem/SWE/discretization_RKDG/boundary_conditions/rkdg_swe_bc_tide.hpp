#ifndef RKDG_SWE_BC_TIDE_HPP
#define RKDG_SWE_BC_TIDE_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace BC {
class Tide {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;

  private:
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    std::vector<DynRowVector<double>> amplitude;
    std::vector<DynRowVector<double>> phase;

    std::vector<DynRowVector<double>> amplitude_gp;
    std::vector<DynRowVector<double>> phase_gp;

  public:
    Tide() = default;
    Tide(const std::vector<TideInput>& tide_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void ComputeFlux(const RKStepper& stepper,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     const HybMatrix<double, SWE::n_variables>& q_in,
                     const HybMatrix<double, SWE::n_variables>& aux_in,
                     HybMatrix<double, SWE::n_variables>& F_hat);
};

Tide::Tide(const std::vector<TideInput>& tide_input) {
    this->frequency    = tide_input[0].frequency;
    this->forcing_fact = tide_input[0].forcing_fact;
    this->equilib_arg  = tide_input[0].equilib_arg;

    uint n_contituents = this->frequency.size();
    uint n_nodes       = tide_input.size();

    this->amplitude.resize(n_contituents);
    this->phase.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude[con].resize(n_nodes);
        this->phase[con].resize(n_nodes);

        for (uint node = 0; node < n_nodes; ++node) {
            this->amplitude[con][node] = tide_input[node].amplitude[con];
            this->phase[con][node]     = tide_input[node].phase[con];
        }
    }
}

template <typename BoundaryType>
void Tide::Initialize(BoundaryType& bound) {
    uint ngp           = bound.data.get_ngp_boundary(bound.bound_id);
    uint n_contituents = this->frequency.size();

    this->q_ex.resize(SWE::n_variables, ngp);

    this->amplitude_gp.resize(n_contituents);
    this->phase_gp.resize(n_contituents);

    for (uint con = 0; con < n_contituents; ++con) {
        this->amplitude_gp[con] = bound.ComputeBoundaryNodalUgp(this->amplitude[con]);
        this->phase_gp[con]     = bound.ComputeBoundaryNodalUgp(this->phase[con]);
    }
}

void Tide::ComputeFlux(const RKStepper& stepper,
                       const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                       const HybMatrix<double, SWE::n_variables>& q_in,
                       const HybMatrix<double, SWE::n_variables>& aux_in,
                       HybMatrix<double, SWE::n_variables>& F_hat) {
    // *** //
    set_constant(this->q_ex, 0.0);

    for (uint con = 0; con < this->frequency.size(); ++con) {
        for (uint gp = 0; gp < columns(q_in); ++gp) {
            row(this->q_ex, SWE::Variables::ze)[gp] += stepper.GetRamp() * this->forcing_fact[con] *
                                                       this->amplitude_gp[con][gp] *
                                                       cos(this->frequency[con] * stepper.GetTimeAtCurrentStage() +
                                                           this->equilib_arg[con] - this->phase_gp[con][gp]);
        }
    }

    row(this->q_ex, SWE::Variables::qx) = row(q_in, SWE::Variables::qx);
    row(this->q_ex, SWE::Variables::qy) = row(q_in, SWE::Variables::qy);

    for (uint gp = 0; gp < columns(q_in); ++gp) {
        column(F_hat, gp) = LLF_flux(
            Global::g, column(q_in, gp), column(this->q_ex, gp), column(aux_in, gp), column(surface_normal, gp));
    }
}
}
}
}

#endif
