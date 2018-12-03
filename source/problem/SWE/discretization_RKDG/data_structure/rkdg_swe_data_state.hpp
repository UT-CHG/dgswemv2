#ifndef RKDG_SWE_DATA_STATE_HPP
#define RKDG_SWE_DATA_STATE_HPP

namespace SWE {
namespace RKDG {
struct StateAccessor {
    std::array<DynRow<double>,SWE::n_variables> q;
    DynRow<double> aux;

    std::array<DynRow<double>,SWE::n_variables> rhs;
    std::array<DynRow<double>,SWE::n_variables> solution;
};

struct StateData {
    using AccessorType = StateAccessor;

    StateData() = default;
    StateData(const uint nelements, const uint ndof) {
        q.fill(DynMatrix<double>(nelements, ndof));
        aux = DynMatrix<double>(nelements, ndof); /* only bath */
        rhs.fill(DynMatrix<double>(nelements, ndof));
        solution.fill(DynMatrix<double>(nelements, ndof));
    }

    std::array<DynMatrix<double>, SWE::n_variables> q;
    DynMatrix<double> aux;

    std::array<DynMatrix<double>, SWE::n_variables> rhs;
    std::array<DynMatrix<double>, SWE::n_variables> solution;

    AccessorType at(const uint index) {
        return AccessorType{make_rows(q,index),
                row(aux,index),
                make_rows(rhs,index),
                make_rows(solution,index)};
    }

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
