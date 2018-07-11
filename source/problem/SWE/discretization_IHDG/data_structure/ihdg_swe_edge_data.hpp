#ifndef IHDG_SWE_EDGE_DATA_HPP
#define IHDG_SWE_EDGE_DATA_HPP

#include "general_definitions.hpp"

#include "ihdg_swe_edge_data_state.hpp"
#include "ihdg_swe_edge_data_internal.hpp"
#include "ihdg_swe_edge_data_global.hpp"
#include "ihdg_swe_edge_data_local.hpp"

namespace SWE {
namespace IHDG {
struct EdgeData {
    EdgeState edge_state;
    EdgeInternal edge_internal;

    EdgeGlobal edge_global;
    EdgeLocal edge_local;

    void initialize() {
        this->edge_state    = EdgeState(this->ndof_global);
        this->edge_internal = EdgeInternal(this->ngp);

        this->edge_global = EdgeGlobal(this->ndof_global, this->ndof_local, this->ngp);
        this->edge_local  = EdgeLocal(this->ndof_global, this->ndof_local);
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
