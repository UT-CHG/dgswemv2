#ifndef LLF_FLUX_HPP
#define LLF_FLUX_HPP

#include "globals.hpp"

#include <array>
#include <cassert>
#include <cmath>

namespace SWE {
//The normal points form the interior side (in) to the exterior side (ex)
inline
void LLF_flux(const double ze_in, const double ze_ex, const double qx_in, const double qx_ex,
               const double qy_in, const double qy_ex, const double bath,
               const std::array<double,2>& normal, double& ze_flux, double& qx_flux, double& qy_flux)
{

  double h_in = ze_in + bath;
  double h_ex = ze_ex + bath;

  assert(h_in >= 0);
  assert(h_ex >= 0);

  //compute normal fluxes
  double qn_in = qx_in * normal[0] + qy_in * normal[1];
  double qn_ex = qx_ex * normal[0] + qy_ex * normal[1];

  double un_in = qn_in/h_in;
  double un_ex = qn_ex/h_ex;

  double max_eigenvalue = std::max( std::abs(un_in) + std::sqrt( Global::g* h_in),
                                    std::abs(un_ex) + std::sqrt( Global::g* h_ex) );


  //compute the potential energy terms
  double pe_in = Global::g * ( std::pow(ze_in,2)/2 + ze_in*bath );
  double pe_ex = Global::g * ( std::pow(ze_ex,2)/2 + ze_ex*bath );

  ze_flux = 0.5*( qn_in + qn_ex + max_eigenvalue* ( ze_in - ze_ex) );
  qx_flux = 0.5*( un_in*qx_in + un_ex*qx_ex + (pe_in + pe_ex)*normal[0] + max_eigenvalue*(qx_in - qx_ex) );
  qy_flux = 0.5*( un_in*qy_in + un_ex*qy_ex + (pe_in + pe_ex)*normal[1] + max_eigenvalue*(qy_in - qy_ex) );
}
}
#endif
