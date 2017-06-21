#ifndef SWE_KERNELS_HPP
#define SWE_KERNELS_HPP

#include "../../general_definitions.h"

#include "LLF_flux.hpp"
#include "source_terms/swe_source_terms.hpp"

namespace SWE {
template<typename ElementType>
void volume_kernel(const Stepper& stepper, ElementType& elt )
{
  const uint rk_stage = stepper.get_stage();

  auto& state = elt->data.state[rk_stage];
  auto& internal = elt->data.internal;
  
  //get state at Gauss points
  elt->ComputeUgp(state.ze, internal.ze_at_gp);
  elt->ComputeUgp(state.qx, internal.qx_at_gp);
  elt->ComputeUgp(state.qy, internal.qy_at_gp); 
  
  //assemble flux
  for ( uint gp = 0; gp < internal.get_n_gp(); ++gp ) {

    internal.water_column_hgt_at_gp[gp] = internal.ze_at_gp[gp] + internal.bath_at_gp[gp];

    internal.ze_flux_at_gp[X][gp] = internal.qx_at_gp[gp];
    internal.ze_flux_at_gp[Y][gp] = internal.qy_at_gp[gp];

    internal.qx_flux_at_gp[X][gp] = std::pow(internal.qx_at_gp[gp],2) / internal.water_column_hgt_at_gp[gp] +
    Global::g * ( 0.5 * std::pow(internal.ze_at_gp[gp],2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
    internal.qx_flux_at_gp[Y][gp] = internal.qx_at_gp[gp]*internal.qy_at_gp[gp]/internal.water_column_hgt_at_gp[gp];

    internal.qy_flux_at_gp[X][gp] = internal.qx_at_gp[gp]*internal.qy_at_gp[gp]/internal.water_column_hgt_at_gp[gp];
    internal.qy_flux_at_gp[Y][gp] = std::pow(internal.qy_at_gp[gp],2)/ internal.water_column_hgt_at_gp[gp] +
    Global::g * ( 0.5 * std::pow(internal.ze_at_gp[gp],2) + internal.ze_at_gp[gp] * internal.bath_at_gp[gp]);
  }

  //skip dof = 0, which is a constant and thus trivially 0
  for ( uint dof = 1; dof < state.get_ndof(); ++dof ) {
    state.rhs_ze[dof] = elt->IntegrationDPhi(X, dof, internal.ze_flux_at_gp[X]) +
      elt->IntegrationDPhi(Y, dof, internal.ze_flux_at_gp[Y]);

    state.rhs_qx[dof] = elt->IntegrationDPhi(X, dof, internal.qx_flux_at_gp[X]) +
      elt->IntegrationDPhi(Y, dof, internal.qx_flux_at_gp[Y]);

    state.rhs_qy[dof] = elt->IntegrationDPhi(X, dof, internal.qy_flux_at_gp[X]) +
      elt->IntegrationDPhi(Y, dof, internal.qy_flux_at_gp[Y]);  }

  SourceTerms::source_kernel(stepper,elt);
};

//compute contributions due to interior eddge
template<typename InterfaceType>
void interface_kernel(const Stepper& stepper, InterfaceType& intface )
{
  const uint rk_stage = stepper.get_stage();

  auto& state_in = intface->data_in.state[rk_stage];
  auto& boundary_in = intface->data_in.boundary;

  auto& state_ex = intface->data_ex.state[rk_stage];
  auto& boundary_ex = intface->data_ex.boundary;

  intface->ComputeUgpIN(state_in.ze, boundary_in.ze_at_gp);
  intface->ComputeUgpIN(state_in.qx, boundary_in.qx_at_gp);
  intface->ComputeUgpIN(state_in.qy, boundary_in.qy_at_gp);

  intface->ComputeUgpEX(state_ex.ze, boundary_ex.ze_at_gp);
  intface->ComputeUgpEX(state_ex.qx, boundary_ex.qx_at_gp);
  intface->ComputeUgpEX(state_ex.qy, boundary_ex.qy_at_gp);


  //assemble numerical fluxes
  for ( uint gp=0; gp < boundary_in.get_n_gp(); ++gp ) {
    LLF_flux( boundary_in.ze_at_gp[gp], boundary_ex.ze_at_gp[gp],
              boundary_in.qx_at_gp[gp], boundary_ex.qx_at_gp[gp],
              boundary_in.qy_at_gp[gp], boundary_ex.qy_at_gp[gp],
              boundary_in.bath_at_gp[gp], intface->surface_normal[gp],
              boundary_in.ze_numerical_flux_at_gp[gp],
              boundary_in.qx_numerical_flux_at_gp[gp],
              boundary_in.qy_numerical_flux_at_gp[gp]
              );
  }

  //now compute contributions to the righthand side
  for ( uint dof = 0; dof < state_in.get_ndof(); ++dof ) {
    state_in.rhs_ze[dof] -= intface->IntegrationPhiIN(dof,
                                              boundary_in.ze_numerical_flux_at_gp);
    state_in.rhs_qx[dof] -= intface->IntegrationPhiIN(dof,
                                              boundary_in.qx_numerical_flux_at_gp);
    state_in.rhs_qy[dof] -= intface->IntegrationPhiIN(dof,
                                              boundary_in.qy_numerical_flux_at_gp);
  }

  for ( uint dof = 0; dof < state_ex.get_ndof(); ++dof ) {
    state_ex.rhs_ze[dof] += intface->IntegrationPhiEX(dof,
                                              boundary_in.ze_numerical_flux_at_gp);
    state_ex.rhs_qx[dof] += intface->IntegrationPhiEX(dof,
                                              boundary_in.qx_numerical_flux_at_gp);
    state_ex.rhs_qy[dof] += intface->IntegrationPhiEX(dof,
                                              boundary_in.qy_numerical_flux_at_gp);
  }
}

template<typename BoundaryType>
void boundary_kernel(const Stepper& stepper, BoundaryType& bound)
{
  const uint rk_stage = stepper.get_stage();

  auto& state = bound->data.state[rk_stage];
  auto& boundary = bound->data.boundary;

  bound->ComputeUgp(state.ze, boundary.ze_at_gp);
  bound->ComputeUgp(state.qx, boundary.qx_at_gp);
  bound->ComputeUgp(state.qy, boundary.qy_at_gp);

  for ( uint gp=0; gp < boundary.get_n_gp(); ++gp ) {
    double ze_ex, qx_ex, qy_ex;

    switch(bound->type){
      case LAND: 
        double n_x, n_y, t_x, t_y, 
        qn_in, qt_in, qn_ex, qt_ex;

        n_x = bound->surface_normal[X][gp];
        n_y = bound->surface_normal[Y][gp];
        t_x = -n_y;
        t_y = n_x;

        qn_in = boundary.qx_at_gp[gp] * n_x + boundary.qy_at_gp[gp] * n_y;
        qt_in = boundary.qx_at_gp[gp] * t_x + boundary.qy_at_gp[gp] * t_y;

        qn_ex = -qn_in;
        qt_ex = qt_in;

        ze_ex = boundary.ze_at_gp[gp];
        qx_ex = qn_ex*n_x + qt_ex*t_x;
        qy_ex = qn_ex*n_y + qt_ex*t_y;
      break;
      case OCEAN:     
        ze_ex = 0; // tidal force      
        qx_ex = boundary.qx_at_gp[gp];
        qy_ex = boundary.qy_at_gp[gp];
      break;
    }

    /*edg.data->get_ex( stepper.get_t_at_curr_stage(),
                      edg.data->ze_at_gp[gp], edg.data->qx_at_gp[gp], edg.data->qy_at_gp[gp],
                      edg.data->bath_at_gp[gp], normal_at_gp, ze_ex, qx_ex, qy_ex);*/

    //    std::cout << "ze_ex: " << ze_ex << "\n";

    LLF_flux( boundary.ze_at_gp[gp], ze_ex,
              boundary.qx_at_gp[gp], qx_ex,
              boundary.qy_at_gp[gp], qy_ex,
              boundary.bath_at_gp[gp], bound->surface_normal[gp],
              boundary.ze_numerical_flux_at_gp[gp],
              boundary.qx_numerical_flux_at_gp[gp],
              boundary.qy_numerical_flux_at_gp[gp]
              );
  }

  //now compute contributions to the righthand side
  for ( uint dof = 0; dof < state.get_ndof(); ++dof ) {
    state.rhs_ze[dof] -= bound->IntegrationPhi(dof,
                                              boundary.ze_numerical_flux_at_gp);
    state.rhs_qx[dof] -= bound->IntegrationPhi(dof,
                                              boundary.qx_numerical_flux_at_gp);
    state.rhs_qy[dof] -= bound->IntegrationPhi(dof,
                                              boundary.qy_numerical_flux_at_gp);
  }
}

template<typename ElementType>
void update_kernel( const Stepper& stepper, ElementType& elt ) {
  const uint rk_stage = stepper.get_stage();
  double dt = stepper.get_dt();
  
  auto& state = elt->data.state;
  auto& curr_state = elt->data.state[rk_stage];
  auto& next_state = elt->data.state[rk_stage+1];

  curr_state.rhs_ze = elt->SolveLSE(curr_state.rhs_ze);
  curr_state.rhs_qx = elt->SolveLSE(curr_state.rhs_qx);
  curr_state.rhs_qy = elt->SolveLSE(curr_state.rhs_qy);

  std::fill(next_state.ze.begin(), next_state.ze.end(), 0);
  std::fill(next_state.qx.begin(), next_state.qx.end(), 0);
  std::fill(next_state.qy.begin(), next_state.qy.end(), 0);

  for ( uint s = 0; s <= rk_stage; ++s ) {
    for ( uint dof = 0; dof < curr_state.get_ndof(); ++dof ) {
    next_state.ze[dof] += stepper.ark[rk_stage][s]*state[s].ze[dof]
      + dt*stepper.brk[rk_stage][s]*state[s].rhs_ze[dof];

    next_state.qx[dof] += stepper.ark[rk_stage][s]*state[s].qx[dof]
      + dt*stepper.brk[rk_stage][s]*state[s].rhs_qx[dof];

    next_state.qy[dof] += stepper.ark[rk_stage][s]*state[s].qy[dof]
      + dt*stepper.brk[rk_stage][s]*state[s].rhs_qy[dof];
    }
  }
};

template<typename ElementType>
void swap_states( const Stepper& stepper, ElementType& elt)
{
  uint n_stages = stepper.get_num_stages();
  auto& state = elt->data.state;

  std::swap(state[0].ze, state[n_stages].ze);
  std::swap(state[0].qx, state[n_stages].qx);
  std::swap(state[0].qy, state[n_stages].qy);

};

}
#endif
