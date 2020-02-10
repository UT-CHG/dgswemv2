#ifndef SWE_DATA_STATE_HPP
#define SWE_DATA_STATE_HPP

namespace SWE {
struct State {
    State() = default;
    State(const uint ndof)
        : q(SWE::n_variables, ndof),
          aux(1, ndof /* only bath */),
          rhs(SWE::n_variables, ndof),
          solution(SWE::n_variables, ndof),
          b_rhs(ndof),
          b_solution(ndof) {}

    HybMatrix<double, SWE::n_variables> q;
    HybMatrix<double, 1> aux;

    HybMatrix<double, SWE::n_variables> rhs;
    HybMatrix<double, SWE::n_variables> solution;

    DynRowVector<double> b_rhs;
    DynRowVector<double> b_solution;

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

#endif
