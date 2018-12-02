#ifndef EHDG_GN_DATA_WET_DRY_HPP
#define EHDG_GN_DATA_WET_DRY_HPP

namespace GN {
namespace EHDG {
struct WetDry : SWE_SIM::WetDry {
    WetDry() = default;
    WetDry(const uint nvrtx) : SWE_SIM::WetDry(nvrtx) {}
};
}
}

#endif
