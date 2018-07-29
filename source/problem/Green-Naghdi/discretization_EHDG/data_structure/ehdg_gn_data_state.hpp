#ifndef EHDG_GN_DATA_STATE_HPP
#define EHDG_GN_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct State {
    State() = default;
    State(const uint ndof)
        : q(GN::n_variables, ndof),
          aux(GN::n_auxiliaries, ndof),
          rhs(GN::n_variables, ndof),
          solution(GN::n_variables, ndof) {}

    HybMatrix<double, GN::n_variables> q;
    HybMatrix<double, GN::n_auxiliaries> aux;

    HybMatrix<double, GN::n_variables> rhs;
    HybMatrix<double, GN::n_variables> solution;
};
}
}

#endif
