#ifndef GN_DATA_WET_DRY_HPP
#define GN_DATA_WET_DRY_HPP

namespace GN {
struct WetDry : SWE::WetDry {
    WetDry() = default;
    WetDry(const uint nvrtx) : SWE::WetDry(nvrtx) {}
};
}

#endif
