#ifndef EHDG_SWE_DATA_WET_DRY_HPP
#define EHDG_SWE_DATA_WET_DRY_HPP

namespace SWE {
namespace EHDG {
struct WetDry {
    WetDry() = default;
    WetDry(const uint nvrtx)
        : q_lin(SWE::n_variables, nvrtx),
          q_at_vrtx(SWE::n_variables, nvrtx),
          bath_at_vrtx(nvrtx),
          h_at_vrtx(nvrtx),
          h_at_vrtx_temp(nvrtx) {}

    bool wet                 = true;
    bool went_completely_dry = false;

    double bath_min;

    HybMatrix<double, SWE::n_variables> q_lin;
    HybMatrix<double, SWE::n_variables> q_at_vrtx;

    std::vector<double> bath_at_vrtx;
    std::vector<double> h_at_vrtx;
    std::vector<double> h_at_vrtx_temp;
};
}
}

#endif
