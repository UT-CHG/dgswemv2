#ifndef SWE_DATA_INTERNAL_HPP
#define SWE_DATA_INTERNAL_HPP

#include "../../../general_definitions.hpp"

namespace SWE {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : ze_flux_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          qx_flux_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          qy_flux_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          tau_s_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          dp_atm_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          dtidal_pot_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          ze_source_term_at_gp(ngp),
          qx_source_term_at_gp(ngp),
          qy_source_term_at_gp(ngp),
          ze_at_gp(ngp),
          qx_at_gp(ngp),
          qy_at_gp(ngp),
          bath_at_gp(ngp),
          h_at_gp(ngp),
          bath_deriv_wrt_x_at_gp(ngp),
          bath_deriv_wrt_y_at_gp(ngp) {}

    std::array<std::vector<double>, 2> ze_flux_at_gp;
    std::array<std::vector<double>, 2> qx_flux_at_gp;
    std::array<std::vector<double>, 2> qy_flux_at_gp;

    std::array<std::vector<double>, 2> tau_s_at_gp;
    std::array<std::vector<double>, 2> dp_atm_at_gp;
    std::array<std::vector<double>, 2> dtidal_pot_at_gp;

    std::vector<double> ze_source_term_at_gp;
    std::vector<double> qx_source_term_at_gp;
    std::vector<double> qy_source_term_at_gp;

    std::vector<double> ze_at_gp;
    std::vector<double> qx_at_gp;
    std::vector<double> qy_at_gp;
    std::vector<double> bath_at_gp;
    std::vector<double> h_at_gp;

    std::vector<double> bath_deriv_wrt_x_at_gp;
    std::vector<double> bath_deriv_wrt_y_at_gp;

#ifdef HAS_HPX
    template<typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template<typename Archive>
void Internal::serialize(Archive& ar, unsigned) {
    ar & ze_flux_at_gp & qx_flux_at_gp & qy_flux_at_gp
        & ze_source_term_at_gp & qx_source_term_at_gp & qy_source_term_at_gp
        & ze_at_gp & qx_at_gp & qy_at_gp & bath_at_gp & h_at_gp
        & bath_deriv_wrt_x_at_gp & bath_deriv_wrt_y_at_gp;
}
#endif
}
#endif
