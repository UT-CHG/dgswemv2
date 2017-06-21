#ifndef SWE_BOUNDARY_CONDITIONS_HPP
#define SWE_BOUNDARY_CONDITIONS_HPP

#include "../../general_definitions.h"

namespace SWE {
  template<typename BoundaryType>
  void land_bc(const Stepper& stepper, const BoundaryType& bound, uint gp, double& ze_ex, double& qx_ex, double& qy_ex) {
  			auto& boundary = bound->data.boundary;
        
        double n_x, n_y, t_x, t_y, qn_in, qt_in, qn_ex, qt_ex;

				n_x = bound->surface_normal[gp][X];
				n_y = bound->surface_normal[gp][Y];
				t_x = -n_y;
				t_y = n_x;

				qn_in = boundary.qx_at_gp[gp] * n_x + boundary.qy_at_gp[gp] * n_y;
				qt_in = boundary.qx_at_gp[gp] * t_x + boundary.qy_at_gp[gp] * t_y;

				qn_ex = -qn_in;
				qt_ex = qt_in;

				ze_ex = boundary.ze_at_gp[gp];
				qx_ex = qn_ex*n_x + qt_ex*t_x;
				qy_ex = qn_ex*n_y + qt_ex*t_y;
  }

  template<typename BoundaryType>
  void tidal_bc(const Stepper& stepper, const BoundaryType& bound, uint gp, double& ze_ex, double& qx_ex, double& qy_ex) {
   			auto& boundary = bound->data.boundary;
		    
        double H_0 = 0.3;
		    if (stepper.get_t_at_curr_stage() < 172800.0) H_0 = 0.3 *  stepper.get_t_at_curr_stage() / 172800.0; //LINEAR RAMPING
		    double H_ocean = H_0*cos(2 * PI *  stepper.get_t_at_curr_stage() / 43200.0); //FOR TESTING M2 TIDAL WAVE WITH PERIOD OF 12HOURS AND AMPLITUDE OF 0.3m

 				ze_ex = H_ocean;      
				qx_ex = boundary.qx_at_gp[gp];
				qy_ex = boundary.qy_at_gp[gp];
  }
}

#endif