#ifndef GN_DATA_STATE_HPP
#define GN_DATA_STATE_HPP

namespace GN {
struct State : SWE::State {
    State() = default;
    State(const uint ndof)
        : SWE::State(ndof),
#ifdef IHDG_SWE
          rhs(SWE::n_variables, ndof),
          solution(SWE::n_variables, ndof),
#endif
          dbath(GN::n_dimensions, ndof),
          dze(GN::n_dimensions, ndof),
          du(GN::n_du_terms, ndof),
          ddu(GN::n_ddu_terms, ndof),
          ddbath(GN::n_ddbath_terms, ndof),
          dddbath(GN::n_dddbath_terms, ndof),
          w1(GN::n_dimensions, ndof) {
    }

    // These are not present in IHDG simulations
#ifdef IHDG_SWE
    HybMatrix<double, SWE::n_variables> rhs;
    HybMatrix<double, SWE::n_variables> solution;
#endif

    HybMatrix<double, GN::n_dimensions> dbath;
    HybMatrix<double, GN::n_dimensions> dze;
    HybMatrix<double, GN::n_du_terms> du;

    HybMatrix<double, GN::n_ddu_terms> ddu;
    HybMatrix<double, GN::n_ddbath_terms> ddbath;

    HybMatrix<double, GN::n_dddbath_terms> dddbath;

    HybMatrix<double, GN::n_dimensions> w1;
};
}

#endif
