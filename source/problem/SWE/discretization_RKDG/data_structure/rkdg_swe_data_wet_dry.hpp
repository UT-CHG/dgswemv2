#ifndef RKDG_SWE_DATA_WET_DRY_HPP
#define RKDG_SWE_DATA_WET_DRY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct WetDry {
    WetDry() = default;
    WetDry(const uint nvrtx)
        : q_lin(SWE::n_variables, nvrtx),
          q_at_vrtx(SWE::n_variables, nvrtx),
          bath_at_vrtx(nvrtx),
          h_at_vrtx(nvrtx),
          h_at_vrtx_temp(nvrtx) {}

    bool wet;
    bool went_completely_dry;

    double bath_min;

    HybMatrix<double, SWE::n_variables> q_lin;
    HybMatrix<double, SWE::n_variables> q_at_vrtx;

    std::vector<double> bath_at_vrtx;
    std::vector<double> h_at_vrtx;
    std::vector<double> h_at_vrtx_temp;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        /*ar  & wet
            & went_completely_dry   
            & bath_min
            & q_lin
            & q_at_vrtx
            & bath_at_vrtx
            & h_at_vrtx
            & h_at_vrtx_temp;*/
        // clang-format on
    }
#endif
};
}
}

#endif
