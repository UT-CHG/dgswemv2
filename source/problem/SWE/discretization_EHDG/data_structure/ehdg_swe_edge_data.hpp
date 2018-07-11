#ifndef EHDG_SWE_EDGE_DATA_HPP
#define EHDG_SWE_EDGE_DATA_HPP

#include "general_definitions.hpp"

#include "ehdg_swe_edge_data_state.hpp"
#include "ehdg_swe_edge_data_internal.hpp"
#include "ehdg_swe_edge_data_global.hpp"

namespace SWE {
namespace EHDG {
struct EdgeData {
    EdgeState edge_state;
    EdgeInternal edge_internal;

    EdgeGlobal edge_global;

    void initialize() {
        this->edge_state    = EdgeState(this->ndof_global);
        this->edge_internal = EdgeInternal(this->ngp);

        this->edge_global = EdgeGlobal(this->ndof_global, this->ngp);
    }

    uint get_ndof_global() { return this->ndof_global; }
    uint get_ngp() { return this->ngp; }

    void set_ndof_global(const uint ndof_global) { this->ndof_global = ndof_global; }
    void set_ndof_local(const uint ndof_local) { this->ndof_local = ndof_local; }
    void set_ngp(const uint ngp) { this->ngp = ngp; }

  private:
    uint ndof_global;
    uint ndof_local;
    uint ngp;
};
}
}

#endif
