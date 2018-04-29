#ifndef SWE_BC_TIDAL_HPP
#define SWE_BC_TIDAL_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace BC {
class Tidal {
  public:
    void ComputeFlux(const Stepper&             stepper,
                     const Array2D<double>&     surface_normal,
                     const std::vector<double>& sp_in,
                     const std::vector<double>& bath_in,
                     const std::vector<double>& ze_in,
                     const std::vector<double>& qx_in,
                     const std::vector<double>& qy_in,
                     std::vector<double>&       ze_numerical_flux,
                     std::vector<double>&       qx_numerical_flux,
                     std::vector<double>&       qy_numerical_flux);

    void GetEX(const Stepper&             stepper,
               const uint                 gp,
               const Array2D<double>&     surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double&                    ze_ex,
               double&                    qx_ex,
               double&                    qy_ex);
};

void Tidal::ComputeFlux(const Stepper&             stepper,
                        const Array2D<double>&     surface_normal,
                        const std::vector<double>& sp_in,
                        const std::vector<double>& bath_in,
                        const std::vector<double>& ze_in,
                        const std::vector<double>& qx_in,
                        const std::vector<double>& qy_in,
                        std::vector<double>&       ze_numerical_flux,
                        std::vector<double>&       qx_numerical_flux,
                        std::vector<double>&       qy_numerical_flux) {
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

void Tidal::GetEX(const Stepper&             stepper,
                  const uint                 gp,
                  const Array2D<double>&     surface_normal,
                  const std::vector<double>& ze_in,
                  const std::vector<double>& qx_in,
                  const std::vector<double>& qy_in,
                  double&                    ze_ex,
                  double&                    qx_ex,
                  double&                    qy_ex) {
    double ze_0   = 0.2;
    double ze_amp = ze_0;

    // ze_amp = ze_0 * tanh(2 * stepper.GetTimeAtCurrentStage() / (0.25 * 86400.0));  // TANH RAMP

    /*if (stepper.GetTimeAtCurrentStage() < 43200.0) {
        ze_amp = ze_0 * stepper.GetTimeAtCurrentStage() / 43200.0;  // LINEAR RAMPING
    } else {
        ze_amp = ze_0;
    }*/

    ze_ex = ze_amp;  //* cos(2 * PI * stepper.GetTimeAtCurrentStage() / 43200.0);  // M2
    qx_ex = qx_in[gp];
    qy_ex = qy_in[gp];
}
}
}

#endif
