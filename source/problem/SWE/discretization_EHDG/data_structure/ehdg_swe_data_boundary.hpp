#ifndef EHDG_SWE_DATA_BOUNDARY_HPP
#define EHDG_SWE_DATA_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : ze_at_gp(ngp),
          qx_at_gp(ngp),
          qy_at_gp(ngp),
          bath_at_gp(ngp),
          h_at_gp(ngp),
          ze_flux_dot_n_at_gp(ngp),
          qx_flux_dot_n_at_gp(ngp),
          qy_flux_dot_n_at_gp(ngp),
          ze_numerical_flux_at_gp(ngp),
          qx_numerical_flux_at_gp(ngp),
          qy_numerical_flux_at_gp(ngp) {}

    std::vector<double> ze_at_gp;
    std::vector<double> qx_at_gp;
    std::vector<double> qy_at_gp;
    std::vector<double> bath_at_gp;
    std::vector<double> h_at_gp;

    std::vector<double> ze_flux_dot_n_at_gp;
    std::vector<double> qx_flux_dot_n_at_gp;
    std::vector<double> qy_flux_dot_n_at_gp;

    std::vector<double> ze_numerical_flux_at_gp;
    std::vector<double> qx_numerical_flux_at_gp;
    std::vector<double> qy_numerical_flux_at_gp;
};
}
}

#endif
