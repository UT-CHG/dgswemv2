#ifndef EHDG_SWE_EDGE_DATA_HPP
#define EHDG_SWE_EDGE_DATA_HPP

#include "general_definitions.hpp"

#include "ehdg_swe_edge_data_state.hpp"
#include "ehdg_swe_edge_data_global.hpp"

namespace SWE {
namespace EHDG {
struct EdgeData {
    EdgeState edge_state;

    EdgeGlobal edge_global;

    EdgeData() = default;
    EdgeData(const uint ndof, const uint ngp)
        : edge_state(EdgeState(ndof, ngp)), edge_global(EdgeGlobal(ndof, ngp)), ndof(ndof), ngp(ngp) {}

    uint get_ndof() { return this->ndof; }
    uint get_ngp() { return this->ngp; }

  private:
    uint ndof;
    uint ngp;
};
}
}

#endif
