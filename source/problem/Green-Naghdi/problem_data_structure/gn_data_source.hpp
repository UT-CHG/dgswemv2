#ifndef GN_DATA_SOURCE_HPP
#define GN_DATA_SOURCE_HPP

namespace GN {
struct Source : SWE::Source {
    Source() = default;
    Source(const uint nnode) : SWE::Source(nnode) {}

    bool dispersive_correction = true;
};
}

#endif