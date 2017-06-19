#ifndef SWE_SOURCE_TERMS_HPP
#define SWE_SOURCE_TERMS_HPP

#include <algorithm>
#include <cmath>
#include <iostream>

#include "globals.hpp"
#include "stepper.hpp"

namespace SWE {
namespace SourceTerms {
template<typename ElementType>
void source_kernel(const Stepper& stepper, ElementType& elt )
{
  //note that dat is a std::unique_ptr type
  auto& dat = elt.data;
  const uint rk_stage = stepper.get_stage();

  //note we assume that the values at gauss points have already been computed
  for ( uint gp = 0; gp < elt.get_num_area_gauss_points(); ++gp ) {
    //compute contribution of hydrostatic pressure
    dat->qx_source_term_at_gp[gp] = Global::g * dat->bath_deriv_wrt_x_at_gp[gp]
      * dat->ze_at_gp[gp];

    dat->qy_source_term_at_gp[gp] = Global::g * dat->bath_deriv_wrt_y_at_gp[gp]
      * dat->ze_at_gp[gp];

    double u_at_gp = dat->qx_at_gp[gp]/dat->water_column_hgt_at_gp[gp];

    double v_at_gp = dat->qy_at_gp[gp]/dat->water_column_hgt_at_gp[gp];

    //compute bottom friction contribution
    double bottom_friction_stress = Global::Cf *
      std::hypot(u_at_gp, v_at_gp)/dat->water_column_hgt_at_gp[gp];

    dat->qx_source_term_at_gp[gp] -= bottom_friction_stress * dat->qx_at_gp[gp];
    dat->qy_source_term_at_gp[gp] -= bottom_friction_stress * dat->qy_at_gp[gp];
  }

  //update the current rk_stage
  for ( uint dof = 0; dof < elt.get_num_dof(); ++dof ) {
    //current we don't have any source terms that affect ze
    //dat->state[rk_stage].rhs_ze[dof] += elt.integrate_area_phi(dof, dat->ze_flux_at_gp);
    dat->state[rk_stage].rhs_qx[dof] += elt.integrate_area_phi(dof, dat->qx_source_term_at_gp);
    dat->state[rk_stage].rhs_qy[dof] += elt.integrate_area_phi(dof, dat->qy_source_term_at_gp);
  }
};
}
}
#endif
