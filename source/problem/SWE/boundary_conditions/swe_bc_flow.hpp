#ifndef SWE_BC_FLOW_HPP
#define SWE_BC_FLOW_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"

namespace SWE {
class Flow {
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
        double qn_0 = -0.75;
        double qn_amp = 0;
        qn_amp = qn_0 * tanh(2 * stepper.GetTimeAtCurrentStage() / (0.5 * 86400.0));  // TANH RAMP

        double n_x, n_y, t_x, t_y, qn_ex, qt_ex;

        n_x = surface_normal[gp][GlobalCoord::x];
        n_y = surface_normal[gp][GlobalCoord::y];
        t_x = -n_y;
        t_y = n_x;

        qn_ex = qn_amp * cos(2 * PI * stepper.GetTimeAtCurrentStage() / 43200.0);  // M2
        qt_ex = 0;

        ze_ex = ze_in[gp];
        qx_ex = qn_ex * n_x + qt_ex * t_x;
        qy_ex = qn_ex * n_y + qt_ex * t_y;
    }
};
}

#endif