#ifndef EHDG_SWE_DATA_INTERNAL_HPP
#define EHDG_SWE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : q_at_gp(ngp),
          aux_at_gp(ngp),
          Fx_at_gp(ngp),
          Fy_at_gp(ngp),
          source_at_gp(ngp),
          tau_s_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          dp_atm_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          dtide_pot_at_gp({std::vector<double>(ngp), std::vector<double>(ngp)}),
          bath_deriv_wrt_x_at_gp(ngp),
          bath_deriv_wrt_y_at_gp(ngp) {}

    std::vector<Vector<double, SWE::n_variables>> q_at_gp;
    std::vector<Vector<double, SWE::n_auxiliaries>> aux_at_gp;

    std::vector<Vector<double, SWE::n_variables>> Fx_at_gp;
    std::vector<Vector<double, SWE::n_variables>> Fy_at_gp;

    std::vector<Vector<double, SWE::n_variables>> source_at_gp;
    std::array<std::vector<double>, 2> tau_s_at_gp;
    std::array<std::vector<double>, 2> dp_atm_at_gp;
    std::array<std::vector<double>, 2> dtide_pot_at_gp;
    std::vector<double> bath_deriv_wrt_x_at_gp;
    std::vector<double> bath_deriv_wrt_y_at_gp;
};
}
}

#endif
