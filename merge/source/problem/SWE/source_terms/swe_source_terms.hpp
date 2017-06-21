#ifndef SWE_SOURCE_TERMS_HPP
#define SWE_SOURCE_TERMS_HPP

#include "../../../general_definitions.h"
#include "../../../stepper.hpp"

namespace SWE {
namespace SourceTerms {
template<typename ElementType>
void source_kernel(const Stepper& stepper, ElementType& elt )
{
  const uint rk_stage = stepper.get_stage();

  auto& state = elt->data.state[rk_stage];
  auto& internal = elt->data.internal;

  //note we assume that the values at gauss points have already been computed
  for ( uint gp = 0; gp < internal.get_n_gp(); ++gp ) {
    //compute contribution of hydrostatic pressure
    internal.qx_source_term_at_gp[gp] = Global::g * internal.bath_deriv_wrt_x_at_gp[gp]
      * internal.ze_at_gp[gp];

    internal.qy_source_term_at_gp[gp] = Global::g * internal.bath_deriv_wrt_y_at_gp[gp]
      * internal.ze_at_gp[gp];

    double u_at_gp = internal.qx_at_gp[gp]/internal.water_column_hgt_at_gp[gp];

    double v_at_gp = internal.qy_at_gp[gp]/internal.water_column_hgt_at_gp[gp];

    //compute bottom friction contribution
    double bottom_friction_stress = Global::Cf *
      std::hypot(u_at_gp, v_at_gp)/internal.water_column_hgt_at_gp[gp];

    internal.qx_source_term_at_gp[gp] -= bottom_friction_stress * internal.qx_at_gp[gp];
    internal.qy_source_term_at_gp[gp] -= bottom_friction_stress * internal.qy_at_gp[gp];
  }

  //update the current rk_stage
  for ( uint dof = 0; dof < state.get_ndof(); ++dof ) {
    //current we don't have any source terms that affect ze
    state.rhs_qx[dof] += elt->IntegrationPhi(dof, internal.qx_source_term_at_gp);
    state.rhs_qy[dof] += elt->IntegrationPhi(dof, internal.qy_source_term_at_gp);
  }
};
}
}
#endif
