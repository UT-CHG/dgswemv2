#ifndef SWE_BC_TIDAL_HPP
#define SWE_BC_TIDAL_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"

namespace SWE {
class Tidal {
  public:
    void GetEX(const Stepper& stepper,
               const uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex) {
        double ze_0 = 0.1;
        double ze_amp = ze_0;

        ze_amp = ze_0 * tanh(2 * stepper.get_t_at_curr_stage() / (0.25 * 86400.0));  // TANH RAMP
        
        /*if (stepper.get_t_at_curr_stage() < 43200.0) {
            ze_amp = ze_0 * stepper.get_t_at_curr_stage() / 43200.0;  // LINEAR RAMPING
        } else {
            ze_amp = ze_0;
        }*/

        ze_ex = ze_amp;  //* cos(2 * PI * stepper.get_t_at_curr_stage() / 43200.0);  // M2
        qx_ex = qx_in[gp];
        qy_ex = qy_in[gp];
    }
};
}

#endif
