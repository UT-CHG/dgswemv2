#ifndef EHDG_GN_DATA_SOURCE_HPP
#define EHDG_GN_DATA_SOURCE_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct Source : SWE::EHDG::Source {
    Source() = default;
    Source(const uint nnode) : SWE::EHDG::Source(nnode) {}
};
}
}

#endif