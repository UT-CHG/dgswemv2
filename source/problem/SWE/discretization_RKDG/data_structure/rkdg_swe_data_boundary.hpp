#ifndef RKDG_SWE_DATA_BOUNDARY_HPP
#define RKDG_SWE_DATA_BOUNDARY_HPP

namespace SWE {
namespace RKDG {
struct BoundaryAccessor {
    BoundaryAccessor() = default;
    BoundaryAccessor( std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp_,
                      std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp_,
                      std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> F_hat_at_gp_)
        :  q_at_gp(q_at_gp_), aux_at_gp(aux_at_gp_), F_hat_at_gp(F_hat_at_gp_)
    {}

    std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp;
    std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp;

    std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> F_hat_at_gp;
};

struct BoundaryData {
    using AccessorType = BoundaryAccessor;

    BoundaryData()=default;
    BoundaryData(const uint n_interfaces, const uint ngp) {
        q_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(n_interfaces, ngp));
        aux_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(n_interfaces, ngp));
        F_hat_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(n_interfaces, ngp));
    }

    AccessorType at(const uint index) {
        return AccessorType( make_rows_as_views(q_at_gp, index),
                             make_rows_as_views(aux_at_gp, index),
                             make_rows_as_views(F_hat_at_gp, index) );
    }

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp;
    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp;

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> F_hat_at_gp;
#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & q_at_gp
            & aux_at_gp
            & F_hat_at_gp;
        // clang-format on
    }
#endif

};
}
}

#endif
