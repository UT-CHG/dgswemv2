#ifndef SWE_BC_TIDE_HPP
#define SWE_BC_TIDE_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/rkdg_simulation/rk_stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace BC {
class Tide {
  private:
    std::vector<double> frequency;
    std::vector<double> forcing_fact;
    std::vector<double> equilib_arg;

    Array2D<double> amplitude;
    Array2D<double> phase;

    Array2D<double> amplitude_gp;
    Array2D<double> phase_gp;

  public:
    Tide() = default;
    Tide(const std::vector<TideInput>& tide_input);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void ComputeFlux(const RKStepper& stepper,
                     const Array2D<double>& surface_normal,
                     const std::vector<double>& sp_in,
                     const std::vector<double>& bath_in,
                     const std::vector<double>& ze_in,
                     const std::vector<double>& qx_in,
                     const std::vector<double>& qy_in,
                     std::vector<double>& ze_numerical_flux,
                     std::vector<double>& qx_numerical_flux,
                     std::vector<double>& qy_numerical_flux);

    void GetEX(const RKStepper& stepper,
               const uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex);
};

Tide::Tide(const std::vector<TideInput>& tide_input) {
    this->frequency    = tide_input[0].frequency;
    this->forcing_fact = tide_input[0].forcing_fact;
    this->equilib_arg  = tide_input[0].equilib_arg;

    uint n_contituents = this->frequency.size();
    uint n_nodes       = tide_input.size();

    this->amplitude.resize(n_contituents);
    this->phase.resize(n_contituents);

    for (uint con = 0; con < n_contituents; con++) {
        this->amplitude[con].resize(n_nodes);
        this->phase[con].resize(n_nodes);

        for (uint node = 0; node < n_nodes; node++) {
            this->amplitude[con][node] = tide_input[node].amplitude[con];
            this->phase[con][node]     = tide_input[node].phase[con];
        }
    }
}

template <typename BoundaryType>
void Tide::Initialize(BoundaryType& bound) {
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

void Tide::ComputeFlux(const RKStepper& stepper,
                       const Array2D<double>& surface_normal,
                       const std::vector<double>& sp_in,
                       const std::vector<double>& bath_in,
                       const std::vector<double>& ze_in,
                       const std::vector<double>& qx_in,
                       const std::vector<double>& qy_in,
                       std::vector<double>& ze_numerical_flux,
                       std::vector<double>& qx_numerical_flux,
                       std::vector<double>& qy_numerical_flux) {
    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < ze_in.size(); ++gp) {
        this->GetEX(stepper, gp, surface_normal, ze_in, qx_in, qy_in, ze_ex, qx_ex, qy_ex);

        LLF_flux(Global::g,
                 ze_in[gp],
                 ze_ex,
                 qx_in[gp],
                 qx_ex,
                 qy_in[gp],
                 qy_ex,
                 bath_in[gp],
                 sp_in[gp],
                 surface_normal[gp],
                 ze_numerical_flux[gp],
                 qx_numerical_flux[gp],
                 qy_numerical_flux[gp]);
    }
}

void Tide::GetEX(const RKStepper& stepper,
                 const uint gp,
                 const Array2D<double>& surface_normal,
                 const std::vector<double>& ze_in,
                 const std::vector<double>& qx_in,
                 const std::vector<double>& qy_in,
                 double& ze_ex,
                 double& qx_ex,
                 double& qy_ex) {
    double ze = 0;

    double frequency;
    double forcing_fact;
    double eq_argument;

    for (uint con = 0; con < this->frequency.size(); con++) {
        frequency    = this->frequency[con];
        forcing_fact = this->forcing_fact[con];
        eq_argument  = this->equilib_arg[con];

        ze += stepper.GetRamp() * forcing_fact * this->amplitude_gp[con][gp] *
              cos(frequency * stepper.GetTimeAtCurrentStage() + eq_argument - this->phase_gp[con][gp]);
    }

    ze_ex = ze;
    qx_ex = qx_in[gp];
    qy_ex = qy_in[gp];
}
}
}

#endif
