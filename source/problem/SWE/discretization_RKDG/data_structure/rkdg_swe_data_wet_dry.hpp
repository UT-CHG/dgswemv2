#ifndef RKDG_SWE_DATA_WET_DRY_HPP
#define RKDG_SWE_DATA_WET_DRY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct WetDry {
    WetDry() = default;
    WetDry(const uint nvrtx)
        : ze_lin(nvrtx),
          qx_lin(nvrtx),
          qy_lin(nvrtx),
          ze_at_vrtx(nvrtx),
          qx_at_vrtx(nvrtx),
          qy_at_vrtx(nvrtx),
          bath_at_vrtx(nvrtx),
          h_at_vrtx(nvrtx),
          h_at_vrtx_temp(nvrtx) {}

    bool wet;
    bool went_completely_dry;

    double bath_min;

    std::vector<double> ze_lin;
    std::vector<double> qx_lin;
    std::vector<double> qy_lin;

    std::vector<double> ze_at_vrtx;
    std::vector<double> qx_at_vrtx;
    std::vector<double> qy_at_vrtx;
    std::vector<double> bath_at_vrtx;
    std::vector<double> h_at_vrtx;
    std::vector<double> h_at_vrtx_temp;

#ifdef HAS_HPX
    template<typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void WetDry::serialize(Archive& ar, unsigned) {
    ar & bath_min & water_volume & ze_lin & qx_lin & qy_lin & ze_at_vrtx & qx_at_vrtx & qy_at_vrtx & bath_at_vrtx
        & h_at_vrtx & h_at_vrtx_temp & wet;
}
#endif
}
}

#endif
