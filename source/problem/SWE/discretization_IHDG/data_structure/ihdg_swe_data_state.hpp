#ifndef IHDG_SWE_DATA_STATE_HPP
#define IHDG_SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct State {
    State() = default;
    State(const uint ndof) : q(SWE::n_variables, ndof), aux(1, ndof /* only bath */) {}

    HybMatrix<double, SWE::n_variables> q;
    HybMatrix<double, 1> aux;
};
}
}

#endif
