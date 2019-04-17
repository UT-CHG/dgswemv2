#ifndef SWE_DATA_STATE_HPP
#define SWE_DATA_STATE_HPP

namespace SWE {
struct StateAccessor {
    std::array<DynView<double, SO::ColumnMajor>,SWE::n_variables> q;
    DynView<double, SO::ColumnMajor> aux;

    std::array<DynView<double, SO::ColumnMajor>,SWE::n_variables> rhs;
    std::array<DynView<double, SO::ColumnMajor>,SWE::n_variables> solution;
};

struct StateData {
    using AccessorType = StateAccessor;

    StateData() = default;
    StateData(const uint nelements, const uint ndof) {
        q.fill(DynMatrix<double>(nelements, ndof));
        aux = DynMatrix<double>(nelements, ndof); /* only bath */
        rhs.fill(DynMatrix<double>(nelements, ndof));
        solution.fill(DynMatrix<double>(nelements, ndof));

        set_constant(aux, 0.0);
        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            set_constant(q[var], 0.0);
            set_constant(rhs[var], 0.0);
            set_constant(solution[var], 0.0);
        }
    }

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> q;
    DynMatrix<double, SO::ColumnMajor> aux;

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> rhs;
    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> solution;

    AccessorType at(const uint index) {
        return AccessorType{make_rows_as_views(q,index),
                row_as_view(aux,index),
                make_rows_as_views(rhs,index),
                make_rows_as_views(solution,index)};
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

#endif
