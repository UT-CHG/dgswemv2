#ifndef RKDG_SWE_DATA_INTERNAL_HPP
#define RKDG_SWE_DATA_INTERNAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct Internal {
    Internal() = default;
    Internal(const uint ngp)
        : q_at_gp(ngp),
          aux_at_gp(ngp),
          Fx_at_gp(ngp),
          Fy_at_gp(ngp),
          source_at_gp(ngp),
          dbath_at_gp(ngp),
          tau_s_at_gp(ngp),
          dp_atm_at_gp(ngp),
          dtide_pot_at_gp(ngp) {}

    std::vector<Vector<double, SWE::n_variables>> q_at_gp;
    std::vector<Vector<double, SWE::n_auxiliaries>> aux_at_gp;

    std::vector<Vector<double, SWE::n_variables>> Fx_at_gp;
    std::vector<Vector<double, SWE::n_variables>> Fy_at_gp;

    std::vector<Vector<double, SWE::n_variables>> source_at_gp;
    std::vector<Vector<double, SWE::n_dimensions>> dbath_at_gp;
    std::vector<Vector<double, SWE::n_dimensions>> tau_s_at_gp;
    std::vector<Vector<double, SWE::n_dimensions>> dp_atm_at_gp;
    std::vector<Vector<double, SWE::n_dimensions>> dtide_pot_at_gp;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        /*ar  & q_at_gp
            & bath_at_gp
            & h_at_gp 
            & Fx_at_gp
            & Fy_at_gp
            & source_at_gp
            & tau_s_at_gp
            & dp_atm_at_gp
            & dtide_pot_at_gp
            & bath_deriv_wrt_x_at_gp
            & bath_deriv_wrt_y_at_gp;*/
        // clang-format on
    }
#endif
};
}
}

#endif
