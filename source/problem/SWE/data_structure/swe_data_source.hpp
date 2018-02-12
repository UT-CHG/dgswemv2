#ifndef SWE_DATA_SOURCE_HPP
#define SWE_DATA_SOURCE_HPP

#include "../../../general_definitions.hpp"

namespace SWE {
struct Source {
    Source() = default;
    Source(const uint nvrtx)
        : tau_s({std::vector<double>(nvrtx), std::vector<double>(nvrtx)}),
          p_atm(nvrtx),
          tidal_pot(nvrtx),
          manning_n(nvrtx) {}

    double coriolis_f = 0.0;

    bool manning = false;
    double g_manning_n_sq = 0.0;

    std::array<std::vector<double>, 2> tau_s;
    std::vector<double> p_atm;
    std::vector<double> tidal_pot;
    std::vector<double> manning_n;
};
}

#endif