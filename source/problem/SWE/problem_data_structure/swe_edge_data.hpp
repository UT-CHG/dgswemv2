#ifndef SWE_EDGE_DATA_HPP
#define SWE_EDGE_DATA_HPP

#include "swe_edge_data_state.hpp"
#include "swe_edge_data_internal.hpp"

namespace SWE {
struct EdgeData {
    SWE::EdgeState edge_state;
    SWE::EdgeInternal edge_internal;

    void initialize() {
        this->edge_state    = SWE::EdgeState(this->ndof);
        this->edge_internal = SWE::EdgeInternal(this->ngp);
    }

    uint get_ndof() { return this->ndof; }
    uint get_ngp() { return this->ngp; }

    void set_ndof(const uint ndof) { this->ndof = ndof; }
    void set_ngp(const uint ngp) { this->ngp = ngp; }

  private:
    uint ndof;
    uint ngp;
};
}

#endif
