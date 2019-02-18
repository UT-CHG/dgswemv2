#ifndef EHDG_GN_DATA_SOURCE_HPP
#define EHDG_GN_DATA_SOURCE_HPP

namespace GN {
namespace EHDG {
struct Source : SWE::Source {
    Source() = default;
    Source(const uint nnode) : SWE::Source(nnode) {}
};
}
}

#endif