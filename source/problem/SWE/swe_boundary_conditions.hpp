#ifndef SWE_BOUNDARY_CONDITIONS_HPP
#define SWE_BOUNDARY_CONDITIONS_HPP

#include "../../simulation/stepper.hpp"

namespace SWE {
class Land {
  public:
    void GetEX(const Stepper& stepper,
               uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex) {
        double n_x, n_y, t_x, t_y, qn_in, qt_in, qn_ex, qt_ex;

        n_x = surface_normal[gp][GlobalCoord::x];
        n_y = surface_normal[gp][GlobalCoord::y];
        t_x = -n_y;
        t_y = n_x;

        qn_in = qx_in[gp] * n_x + qy_in[gp] * n_y;
        qt_in = qx_in[gp] * t_x + qy_in[gp] * t_y;

        qn_ex = -qn_in;
        qt_ex = qt_in;

        ze_ex = ze_in[gp];
        qx_ex = qn_ex * n_x + qt_ex * t_x;
        qy_ex = qn_ex * n_y + qt_ex * t_y;
    }
};

class Tidal {
  public:
    void GetEX(const Stepper& stepper,
               uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex) {
        double H_0 = 0.2763;

        // double H_ocean = H_0 * tanh(2 * stepper.get_t_at_curr_stage() / (0.25 * 86400.0)); //FOR TESTING TANH RAMP
        if (stepper.get_t_at_curr_stage() < 172800.0)
            H_0 = H_0 * stepper.get_t_at_curr_stage() / 172800.0;  // LINEAR RAMPING
        double H_ocean = H_0 * cos(2 * PI * stepper.get_t_at_curr_stage() /
                                   43200.0);  // FOR TESTING M2 TIDAL WAVE WITH PERIOD OF 12HOURS

        ze_ex = H_ocean;
        qx_ex = qx_in[gp];
        qy_ex = qy_in[gp];
    }
};
}

#endif