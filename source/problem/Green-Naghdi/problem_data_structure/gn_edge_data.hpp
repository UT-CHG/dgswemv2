#ifndef GN_EDGE_DATA_HPP
#define GN_EDGE_DATA_HPP

#include "gn_edge_data_state.hpp"
#include "gn_edge_data_internal.hpp"

namespace GN {
struct EdgeData {
    GN::EdgeState edge_state;
    GN::EdgeInternal edge_internal;

    void initialize() {
        this->edge_state    = GN::EdgeState(this->ndof);
        this->edge_internal = GN::EdgeInternal(this->ngp);
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
