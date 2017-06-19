#ifndef SWE_KERNELS_HPP
#define SWE_KERNELS_HPP

#include <algorithm>
#include <cmath>
#include <iostream>

#include "globals.hpp"
#include "stepper.hpp"

#include "LLF_flux.hpp"
#include "source_terms/swe_source_terms.hpp"

namespace SWE {
template<typename ElementType>
void volume_kernel(const Stepper& stepper, ElementType& elt )
{
  using Global::g;
  //note that dat is a std::unique_ptr type
  auto& dat = elt.data;

  const uint rk_stage = stepper.get_stage();

  //get state at Gauss points
  elt.fill_area_gauss_points(dat->state[rk_stage].ze, dat->ze_at_gp);
  elt.fill_area_gauss_points(dat->state[rk_stage].qx, dat->qx_at_gp);
  elt.fill_area_gauss_points(dat->state[rk_stage].qy, dat->qy_at_gp);

  //assemble flux
  for ( uint i = 0; i < elt.get_num_area_gauss_points(); ++i ) {
    dat->water_column_hgt_at_gp[i] = dat->ze_at_gp[i] + dat->bath_at_gp[i];

    dat->ze_flux_at_gp[i][0] = dat->qx_at_gp[i];
    dat->ze_flux_at_gp[i][1] = dat->qy_at_gp[i];

    dat->qx_flux_at_gp[i][0] = std::pow(dat->qx_at_gp[i],2)/ dat->water_column_hgt_at_gp[i] +
      g * ( 0.5 * std::pow(dat->ze_at_gp[i],2) + dat->ze_at_gp[i] * dat->bath_at_gp[i]);
    dat->qx_flux_at_gp[i][1] = dat->qx_at_gp[i]*dat->qy_at_gp[i]/dat->water_column_hgt_at_gp[i];

    dat->qy_flux_at_gp[i][0] = dat->qx_at_gp[i]*dat->qy_at_gp[i]/dat->water_column_hgt_at_gp[i];
    dat->qy_flux_at_gp[i][1] = std::pow(dat->qy_at_gp[i],2)/ dat->water_column_hgt_at_gp[i] +
      g * ( 0.5 * std::pow(dat->ze_at_gp[i],2) + dat->ze_at_gp[i] * dat->bath_at_gp[i]);
  }

  //skip dof = 0, which is a constant and thus trivially 0
  for ( uint dof = 1; dof < elt.get_num_dof(); ++dof ) {
    dat->state[rk_stage].rhs_ze[dof] = elt.integrate_area_Dphi(dof, dat->ze_flux_at_gp);
    dat->state[rk_stage].rhs_qx[dof] = elt.integrate_area_Dphi(dof, dat->qx_flux_at_gp);
    dat->state[rk_stage].rhs_qy[dof] = elt.integrate_area_Dphi(dof, dat->qy_flux_at_gp);
  }


  SourceTerms::source_kernel(stepper,elt);
};

//compute contributions due to interior eddge
template<typename EdgeType>
void edge_kernel(const Stepper& stepper, EdgeType& edg )
{
  auto& data_in = edg._elts.first->data;
  auto& data_ex = edg._elts.second->data;

  const uint rk_stage = stepper.get_stage();

  edg.fill_surface_gauss_points(0,data_in->state[rk_stage].ze, edg.data->ze_at_gp[0]);
  edg.fill_surface_gauss_points(1,data_ex->state[rk_stage].ze, edg.data->ze_at_gp[1]);

  edg.fill_surface_gauss_points(0,data_in->state[rk_stage].qx, edg.data->qx_at_gp[0]);
  edg.fill_surface_gauss_points(1,data_ex->state[rk_stage].qx, edg.data->qx_at_gp[1]);

  edg.fill_surface_gauss_points(0,data_in->state[rk_stage].qy, edg.data->qy_at_gp[0]);
  edg.fill_surface_gauss_points(1,data_ex->state[rk_stage].qy, edg.data->qy_at_gp[1]);

  //assemble numerical fluxes
  for ( uint gp=0; gp < edg.get_num_surface_gauss_points(); ++gp ) {

    const std::array<double,2> normal_at_gp = edg.get_normal_at_gp(gp);

    LLF_flux( edg.data->ze_at_gp[0][gp], edg.data->ze_at_gp[1][gp],
              edg.data->qx_at_gp[0][gp], edg.data->qx_at_gp[1][gp],
              edg.data->qy_at_gp[0][gp], edg.data->qy_at_gp[1][gp],
              edg.data->bath_at_gp[gp], normal_at_gp,
              edg.data->ze_numerical_flux_at_gp[gp],
              edg.data->qx_numerical_flux_at_gp[gp],
              edg.data->qy_numerical_flux_at_gp[gp]
              );
  }

  //now compute contributions to the righthand side
  for ( uint dof = 0; dof < edg._elts.first->get_num_dof(); ++dof ) {
    data_in->state[rk_stage].rhs_ze[dof] -= edg.integrate_surface_phi(0, dof,
                                              edg.data->ze_numerical_flux_at_gp);
    data_in->state[rk_stage].rhs_qx[dof] -= edg.integrate_surface_phi(0, dof,
                                              edg.data->qx_numerical_flux_at_gp);
    data_in->state[rk_stage].rhs_qy[dof] -= edg.integrate_surface_phi(0, dof,
                                              edg.data->qy_numerical_flux_at_gp);
  }

  for ( uint dof = 0; dof < edg._elts.second->get_num_dof(); ++dof ) {
    data_ex->state[rk_stage].rhs_ze[dof] += edg.integrate_surface_phi(1, dof,
                                              edg.data->ze_numerical_flux_at_gp);
    data_ex->state[rk_stage].rhs_qx[dof] += edg.integrate_surface_phi(1, dof,
                                              edg.data->qx_numerical_flux_at_gp);
    data_ex->state[rk_stage].rhs_qy[dof] += edg.integrate_surface_phi(1, dof,
                                              edg.data->qy_numerical_flux_at_gp);
  }
}

template<typename EdgeType>
void boundary_kernel(const Stepper& stepper, EdgeType& edg)
{
  auto& data_in = edg._elt->data;
  const uint rk_stage = stepper.get_stage();

  edg.fill_surface_gauss_points(0,data_in->state[rk_stage].ze, edg.data->ze_at_gp);
  edg.fill_surface_gauss_points(0,data_in->state[rk_stage].qx, edg.data->qx_at_gp);
  edg.fill_surface_gauss_points(0,data_in->state[rk_stage].qy, edg.data->qy_at_gp);

  //assemble numerical fluxes
  for ( uint gp=0; gp < edg.get_num_surface_gauss_points(); ++gp ) {

    const std::array<double,2> normal_at_gp = edg.get_normal_at_gp(gp);

    double ze_ex, qx_ex, qy_ex;

    edg.data->get_ex( stepper.get_t_at_curr_stage(),
                      edg.data->ze_at_gp[gp], edg.data->qx_at_gp[gp], edg.data->qy_at_gp[gp],
                      edg.data->bath_at_gp[gp], normal_at_gp, ze_ex, qx_ex, qy_ex);

    //    std::cout << "ze_ex: " << ze_ex << "\n";

    LLF_flux( edg.data->ze_at_gp[gp], ze_ex,
              edg.data->qx_at_gp[gp], qx_ex,
              edg.data->qy_at_gp[gp], qy_ex,
              edg.data->bath_at_gp[gp], normal_at_gp,
              edg.data->ze_numerical_flux_at_gp[gp],
              edg.data->qx_numerical_flux_at_gp[gp],
              edg.data->qy_numerical_flux_at_gp[gp]
              );
  }

  //now compute contributions to the righthand side
  for ( uint dof = 0; dof < edg._elt->get_num_dof(); ++dof ) {
    data_in->state[rk_stage].rhs_ze[dof] -= edg.integrate_surface_phi(0,dof,
                                              edg.data->ze_numerical_flux_at_gp);

    // std::cout << "Dof: " << dof << " rhsZe: " << data_in->state[rk_stage].rhs_ze[dof] << "\n";

    data_in->state[rk_stage].rhs_qx[dof] -= edg.integrate_surface_phi(0,dof,
                                              edg.data->qx_numerical_flux_at_gp);
    data_in->state[rk_stage].rhs_qy[dof] -= edg.integrate_surface_phi(0,dof,
                                              edg.data->qy_numerical_flux_at_gp);
  }
}

template<typename ElementType>
void update_kernel( const Stepper& stepper, ElementType& elt ) {

  const uint rk_stage = stepper.get_stage();
  auto& state = elt.data->state;
  double dt = stepper.get_dt();

  elt.apply_Minv(state[rk_stage].rhs_ze);
  elt.apply_Minv(state[rk_stage].rhs_qx);
  elt.apply_Minv(state[rk_stage].rhs_qy);

  std::fill(state[rk_stage+1].ze.begin(), state[rk_stage+1].ze.end(), 0);
  std::fill(state[rk_stage+1].qx.begin(), state[rk_stage+1].qx.end(), 0);
  std::fill(state[rk_stage+1].qy.begin(), state[rk_stage+1].qy.end(), 0);

  for ( uint s = 0; s <= rk_stage; ++s ) {
    for ( uint dof = 0; dof < elt.get_num_dof(); ++dof ) {
    state[rk_stage+1].ze[dof] += stepper.ark[rk_stage][s]*state[s].ze[dof]
      + dt*stepper.brk[rk_stage][s]*state[s].rhs_ze[dof];

    state[rk_stage+1].qx[dof] += stepper.ark[rk_stage][s]*state[s].qx[dof]
      + dt*stepper.brk[rk_stage][s]*state[s].rhs_qx[dof];

    state[rk_stage+1].qy[dof] += stepper.ark[rk_stage][s]*state[s].qy[dof]
      + dt*stepper.brk[rk_stage][s]*state[s].rhs_qy[dof];
    }
  }
};

template<typename ElementType>
void swap_states( const Stepper& stepper, ElementType& elt)
{
  uint n_stages = stepper.get_num_stages();
  auto& dat = elt.data;

  std::swap(dat->state[0].ze, dat->state[n_stages].ze);
  std::swap(dat->state[0].qx, dat->state[n_stages].qx);
  std::swap(dat->state[0].qy, dat->state[n_stages].qy);

};

}
#endif
