#ifndef GN_DATA_SLOPE_LIMIT_HPP
#define GN_DATA_SLOPE_LIMIT_HPP

namespace GN {
struct SlopeLimit : SWE::SlopeLimit {
    SlopeLimit() = default;
    SlopeLimit(const uint nvrtx, const uint nbound) : SWE::SlopeLimit(nvrtx, nbound) {}
};
}

#endif
