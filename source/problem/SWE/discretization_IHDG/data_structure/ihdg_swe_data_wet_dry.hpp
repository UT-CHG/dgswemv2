#ifndef IHDG_SWE_DATA_WET_DRY_HPP
#define IHDG_SWE_DATA_WET_DRY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct WetDry {
    bool wet                 = true;
    bool went_completely_dry = false;
};
}
}

#endif
