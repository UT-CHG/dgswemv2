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

    HybMatrix<double, SWE::n_variables> q;
    HybMatrix<double, 1> aux;

    HybMatrix<double, SWE::n_variables> rhs;
    HybMatrix<double, SWE::n_variables> solution;
};
}
}

#endif
