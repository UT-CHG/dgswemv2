#ifndef IHDG_SWE_EDGE_DATA_LOCAL_HPP
#define IHDG_SWE_EDGE_DATA_LOCAL_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct EdgeLocal {
    EdgeLocal() = default;
    EdgeLocal(const uint ndof_global, const uint ndof_local)
        : delta_hat_matrix(SWE::n_variables * ndof_local, SWE::n_variables * ndof_global) {}

    DMatrix<double> delta_hat_matrix;
};
}
}

#endif
