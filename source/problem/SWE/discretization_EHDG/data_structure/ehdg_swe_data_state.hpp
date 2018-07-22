#ifndef EHDG_SWE_DATA_STATE_HPP
#define EHDG_SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace EHDG {
struct State {
    State() = default;
    State(const uint ndof)
        : q(SWE::n_variables, ndof),
          aux(1, ndof /* only bath */),
          rhs(SWE::n_variables, ndof),
          solution(SWE::n_variables, ndof) {}

    DynMatrix<double> q;
    DynMatrix<double> aux;

    DynMatrix<double> rhs;
    DynMatrix<double> solution;
};
}
}

#endif
