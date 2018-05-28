#ifndef SWE_BC_TIDAL_HPP
#define SWE_BC_TIDAL_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace BC {
class Tidal {
  private:
    Array2D<double> tidal_constituents;  // 0 - freq, 1 - force_fact, 2 - eq_arg
    Array2D<double> amplitude;
    Array2D<double> phase;

    Array2D<double> amplitude_gp;
    Array2D<double> phase_gp;

  public:
    Tidal() = default;
    Tidal(const Array2D<double>& tidal_constituents, const Array2D<double>& amplitude, const Array2D<double>& phase);

    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void ComputeFlux(const Stepper& stepper,
                     const Array2D<double>& surface_normal,
                     const std::vector<double>& sp_in,
                     const std::vector<double>& bath_in,
                     const std::vector<double>& ze_in,
                     const std::vector<double>& qx_in,
                     const std::vector<double>& qy_in,
                     std::vector<double>& ze_numerical_flux,
                     std::vector<double>& qx_numerical_flux,
                     std::vector<double>& qy_numerical_flux);

    void GetEX(const Stepper& stepper,
               const uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex);
};

Tidal::Tidal(const Array2D<double>& tidal_constituents, const Array2D<double>& amplitude, const Array2D<double>& phase)
    : tidal_constituents(tidal_constituents), amplitude(amplitude), phase(phase) {}

template <typename BoundaryType>
void Tidal::Initialize(BoundaryType& bound) {
    this->amplitude_gp.resize(this->tidal_constituents.size());
    this->phase_gp.resize(this->tidal_constituents.size());

    for (uint con = 0; con < this->tidal_constituents.size(); con++) {
        this->amplitude_gp[con].resize(bound.data.get_ngp_boundary(bound.bound_id));
        this->phase_gp[con].resize(bound.data.get_ngp_boundary(bound.bound_id));

        bound.ComputeBoundaryNodalUgp(this->amplitude[con], this->amplitude_gp[con]);
        bound.ComputeBoundaryNodalUgp(this->phase[con], this->phase_gp[con]);
    }
}

void Tidal::ComputeFlux(const Stepper& stepper,
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

void Tidal::GetEX(const Stepper& stepper,
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

    for (uint con = 0; con < this->tidal_constituents.size(); con++) {
        frequency    = this->tidal_constituents[con][0];
        forcing_fact = this->tidal_constituents[con][1];
        eq_argument  = this->tidal_constituents[con][2];

        ze += forcing_fact * this->amplitude_gp[con][gp] *
              cos(frequency * stepper.GetTimeAtCurrentStage() + eq_argument - this->phase_gp[con][gp]);
    }

    ze_ex = ze;
    qx_ex = qx_in[gp];
    qy_ex = qy_in[gp];
}
}
}

#endif
