#ifndef RKDG_SWE_DATA_INTERNAL_HPP
#define RKDG_SWE_DATA_INTERNAL_HPP

namespace SWE {
namespace RKDG {
struct InternalAccessor {
    InternalAccessor() = default;
    InternalAccessor(std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp_,
                     std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp_)
        : q_at_gp(std::move(q_at_gp_)), aux_at_gp(std::move(aux_at_gp_)) {

        uint ngp     = q_at_gp[0].size();
        Fx_at_gp     = HybMatrix<double, SWE::n_variables>(SWE::n_variables, ngp);
        Fy_at_gp     = HybMatrix<double, SWE::n_variables>(SWE::n_variables, ngp);
        source_at_gp = HybMatrix<double, SWE::n_variables>(SWE::n_variables, ngp);
        dbath_at_gp  = HybMatrix<double, SWE::n_dimensions>(SWE::n_dimensions, ngp);
        tau_s_at_gp  = HybMatrix<double, SWE::n_dimensions>(SWE::n_dimensions, ngp);
        dp_atm_at_gp = HybMatrix<double, SWE::n_dimensions>(SWE::n_dimensions, ngp);
        dtide_pot_at_gp = HybMatrix<double, SWE::n_dimensions>(SWE::n_dimensions, ngp);

    }

    std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp;
    std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp;

    HybMatrix<double, SWE::n_variables> Fx_at_gp;
    HybMatrix<double, SWE::n_variables> Fy_at_gp;

    HybMatrix<double, SWE::n_variables> source_at_gp;
    HybMatrix<double, SWE::n_dimensions> dbath_at_gp;
    HybMatrix<double, SWE::n_dimensions> tau_s_at_gp;
    HybMatrix<double, SWE::n_dimensions> dp_atm_at_gp;
    HybMatrix<double, SWE::n_dimensions> dtide_pot_at_gp;
};

struct InternalData {
    using AccessorType = InternalAccessor;

    InternalData() = default;
    InternalData(const uint nelements, const uint ngp) {
        q_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(nelements,ngp));
        aux_at_gp.fill(DynMatrix<double, SO::ColumnMajor>(nelements,ngp));
        F_at_gp = DynMatrix<double, SO::ColumnMajor>(nelements,ngp);
    }

    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp;
    std::array<DynMatrix<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp;
    StatVector<DynMatrix<double, SO::ColumnMajor>, SWE::n_dimensions> F_at_gp;

    AccessorType at(const uint index) {
        return AccessorType( make_rows_as_views(q_at_gp, index),
                             make_rows_as_views(aux_at_gp, index));
    }

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & q_at_gp
            & aux_at_gp
            & Fx_at_gp
            & Fy_at_gp
            & source_at_gp
            & dbath_at_gp
            & tau_s_at_gp
            & dp_atm_at_gp
            & dtide_pot_at_gp;
        // clang-format on
    }
#endif
};
}
}

#endif
