#ifndef IHDG_SWE_DATA_STATE_HPP
#define IHDG_SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace IHDG {
struct State {
    State() = default;
    State(const uint ndof) : q(SWE::n_variables, ndof), aux(1, ndof /* only bath */) {}

    DynMatrix<double> q;
    DynMatrix<double> aux;
};
}
}

#endif
