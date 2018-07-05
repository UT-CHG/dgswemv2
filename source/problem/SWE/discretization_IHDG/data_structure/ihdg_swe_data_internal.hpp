#ifndef IHDG_SWE_DATA_INTERNAL_HPP
#define IHDG_SWE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : ze_at_gp(ngp),
          qx_at_gp(ngp),
          qy_at_gp(ngp),
          bath_at_gp(ngp),
          h_at_gp(ngp),
          ze_back_at_gp(ngp),
          qx_back_at_gp(ngp),
          qy_back_at_gp(ngp),
          h_back_at_gp(ngp),
          tau_s_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          dp_atm_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          dtide_pot_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          ze_source_term_at_gp(ngp),
          qx_source_term_at_gp(ngp),
          qy_source_term_at_gp(ngp),
          bath_deriv_wrt_x_at_gp(ngp),
          bath_deriv_wrt_y_at_gp(ngp) {}

    std::vector<double> ze_at_gp;
    std::vector<double> qx_at_gp;
    std::vector<double> qy_at_gp;
    std::vector<double> bath_at_gp;
    std::vector<double> h_at_gp;

    std::vector<double> ze_back_at_gp;
    std::vector<double> qx_back_at_gp;
    std::vector<double> qy_back_at_gp;
    std::vector<double> h_back_at_gp;

    std::array<std::vector<double>, 2> tau_s_at_gp;
    std::array<std::vector<double>, 2> dp_atm_at_gp;
    std::array<std::vector<double>, 2> dtide_pot_at_gp;

    std::vector<double> ze_source_term_at_gp;
    std::vector<double> qx_source_term_at_gp;
    std::vector<double> qy_source_term_at_gp;

    std::vector<double> bath_deriv_wrt_x_at_gp;
    std::vector<double> bath_deriv_wrt_y_at_gp;
};
}
}

#endif
