#ifndef RKDG_SWE_DATA_STATE_HPP
#define RKDG_SWE_DATA_STATE_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
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

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        /*ar  & q
            & bath
            & rhs
            & solution;*/
        // clang-format on
    }
#endif
};
}
}

#endif
