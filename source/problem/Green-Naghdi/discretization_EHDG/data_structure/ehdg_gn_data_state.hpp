#ifndef EHDG_GN_DATA_STATE_HPP
#define EHDG_GN_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
struct State {
    State() = default;
    State(const uint ndof)
        : q(GN::n_variables, ndof),
          aux(1, ndof /* only bath */),
          rhs(GN::n_variables, ndof),
          solution(GN::n_variables, ndof),
          dze(GN::n_dimensions, ndof),
          du(GN::n_du_terms, ndof),
          ddu(GN::n_ddu_terms, ndof),
          w1(GN::n_dimensions, ndof) {}

    /* swe containers */

    HybMatrix<double, GN::n_variables> q;
    HybMatrix<double, 1> aux;

    HybMatrix<double, GN::n_variables> rhs;
    HybMatrix<double, GN::n_variables> solution;

    /* dispersive correction containers */

    HybMatrix<double, GN::n_dimensions> dze;
    HybMatrix<double, GN::n_du_terms> du;
    HybMatrix<double, GN::n_ddu_terms> ddu;

    HybMatrix<double, GN::n_dimensions> w1;
};
}
}

#endif
