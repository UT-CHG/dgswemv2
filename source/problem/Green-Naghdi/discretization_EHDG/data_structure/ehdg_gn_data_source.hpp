#ifndef EHDG_GN_DATA_SOURCE_HPP
#define EHDG_GN_DATA_SOURCE_HPP

namespace GN {
namespace EHDG {
struct Source : SWE_SIM::Source {
    Source() = default;
    Source(const uint nnode) : SWE_SIM::Source(nnode) {}
};
}
}

#endif