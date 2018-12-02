#ifndef EHDG_GN_DATA_SLOPE_LIMIT_HPP
#define EHDG_GN_DATA_SLOPE_LIMIT_HPP

namespace GN {
namespace EHDG {
struct SlopeLimit : SWE_SIM::SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nvrtx, const uint nbound) : SWE_SIM::SlopeLimit(nvrtx, nbound) {}
};
}
}

#endif
