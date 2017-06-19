#ifndef SWE_DATA_HPP
#define SWE_DATA_HPP

#include "../../general_definitions.h"

namespace SWE {
struct State {
  State(int ndof)
    : ze(ndof), qx(ndof), qy(ndof),
      rhs_ze(ndof), rhs_qx(ndof), rhs_qy(ndof)
  {}

  std::vector<double> ze;
  std::vector<double> qx;
  std::vector<double> qy;
  std::vector<double> rhs_ze;
  std::vector<double> rhs_qx;
  std::vector<double> rhs_qy;
};

struct Boundary
{
  Boundary(int n_gp)
    : ze_at_gp({ std::vector<double>(n_gp), std::vector<double>(n_gp) }),
      qx_at_gp({ std::vector<double>(n_gp), std::vector<double>(n_gp) }),
      qy_at_gp({ std::vector<double>(n_gp), std::vector<double>(n_gp) }),
      bath_at_gp(n_gp),
      ze_numerical_flux_at_gp(n_gp),
      qx_numerical_flux_at_gp(n_gp),
      qy_numerical_flux_at_gp(n_gp)
  {}

  std::array<std::vector<double>,2> ze_at_gp;
  std::array<std::vector<double>,2> qx_at_gp;
  std::array<std::vector<double>,2> qy_at_gp;

  std::vector<double> bath_at_gp;

  std::vector<double> ze_numerical_flux_at_gp;
  std::vector<double> qx_numerical_flux_at_gp;
  std::vector<double> qy_numerical_flux_at_gp;
};

struct Internal {
  Internal(int n_gp)
    : ze_flux_at_gp(n_gp),
      qx_flux_at_gp(n_gp),
      qy_flux_at_gp(n_gp),
      ze_source_term_at_gp(n_gp),
      qx_source_term_at_gp(n_gp),
      qy_source_term_at_gp(n_gp),
      ze_at_gp(n_gp),
      qx_at_gp(n_gp),
      qy_at_gp(n_gp),
      bath_at_gp(n_gp),
      bath_deriv_wrt_x_at_gp(n_gp),
      bath_deriv_wrt_y_at_gp(n_gp),
      water_column_hgt_at_gp(n_gp)
  {}

  std::vector<std::array<double,2> > ze_flux_at_gp;
  std::vector<std::array<double,2> > qx_flux_at_gp;
  std::vector<std::array<double,2> > qy_flux_at_gp;

  std::vector<double> ze_source_term_at_gp;
  std::vector<double> qx_source_term_at_gp;
  std::vector<double> qy_source_term_at_gp;

  std::vector<double> ze_at_gp;
  std::vector<double> qx_at_gp;
  std::vector<double> qy_at_gp;

  std::vector<double> bath_at_gp;
  std::vector<double> bath_deriv_wrt_x_at_gp;
  std::vector<double> bath_deriv_wrt_y_at_gp;

  std::vector<double> water_column_hgt_at_gp;
};

struct Data {
  Data()
   : state_data(0, State(0)), internal_data(0), boundary_data(0) 
  {}

  std::vector<State> state_data;
  Internal internal_data;
  Boundary boundary_data;
};
}

#endif