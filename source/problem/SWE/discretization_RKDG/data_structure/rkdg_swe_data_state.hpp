#ifndef RKDG_SWE_DATA_STATE_HPP
#define RKDG_SWE_DATA_STATE_HPP

namespace SWE {
namespace RKDG {
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

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & q
            & aux
            & rhs
            & solution;
        // clang-format on
    }
#endif
};
}
}

#endif
