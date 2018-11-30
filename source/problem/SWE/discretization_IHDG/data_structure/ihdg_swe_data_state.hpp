#ifndef IHDG_SWE_DATA_STATE_HPP
#define IHDG_SWE_DATA_STATE_HPP

namespace SWE {
namespace IHDG {
struct State {
    State() = default;
    State(const uint ndof) : q(SWE::n_variables, ndof), aux(1, ndof /* only bath */), rhs(SWE::n_variables, ndof) {}

    HybMatrix<double, SWE::n_variables> q;
    HybMatrix<double, 1> aux;

    HybMatrix<double, SWE::n_variables> rhs;
};
}
}

#endif
