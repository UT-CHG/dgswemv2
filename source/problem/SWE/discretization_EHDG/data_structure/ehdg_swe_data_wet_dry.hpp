#ifndef EHDG_SWE_DATA_WET_DRY_HPP
#define EHDG_SWE_DATA_WET_DRY_HPP

namespace SWE {
namespace EHDG {
struct WetDry {
    bool wet                 = true;
    bool went_completely_dry = false;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & wet
            & went_completely_dry;
        // clang-format on
    }
#endif
};
}
}

#endif
