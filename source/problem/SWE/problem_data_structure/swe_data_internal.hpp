#ifndef SWE_DATA_INTERNAL_HPP
#define SWE_DATA_INTERNAL_HPP

namespace SWE {
struct InternalAccessor {
    InternalAccessor() = default;
    InternalAccessor(std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp_,
                     std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp_)
        : q_at_gp(std::move(q_at_gp_)), aux_at_gp(std::move(aux_at_gp_)) {

        uint ngp     = q_at_gp[0].size();
        Fx_at_gp     = DynMatrix<double>(SWE::n_variables, ngp);
        Fy_at_gp     = DynMatrix<double>(SWE::n_variables, ngp);
        source_at_gp = DynMatrix<double>(SWE::n_variables, ngp);
        dbath_at_gp  = DynMatrix<double>(SWE::n_dimensions, ngp);
        tau_s_at_gp  = DynMatrix<double>(SWE::n_dimensions, ngp);
        dp_atm_at_gp = DynMatrix<double>(SWE::n_dimensions, ngp);
        dtide_pot_at_gp = DynMatrix<double>(SWE::n_dimensions, ngp);

        q_prev_at_gp   = HybMatrix<double, SWE::n_variables>(SWE::n_variables, ngp);
        del_q_DT_at_gp = HybMatrix<double, SWE::n_variables>(SWE::n_variables, ngp);
        kronecker_DT_at_gp = HybMatrix<double, SWE::n_variables * SWE::n_variables>(SWE::n_variables * SWE::n_variables, ngp);
        dFx_dq_at_gp   = HybMatrix<double, SWE::n_variables * SWE::n_variables>(SWE::n_variables * SWE::n_variables, ngp);
        dFy_dq_at_gp   = HybMatrix<double, SWE::n_variables * SWE::n_variables>(SWE::n_variables * SWE::n_variables, ngp);
        dsource_dq_at_gp = HybMatrix<double, SWE::n_variables * SWE::n_variables>(SWE::n_variables * SWE::n_variables, ngp);

    }

    std::array<DynView<double, SO::ColumnMajor>, SWE::n_variables> q_at_gp;
    std::array<DynView<double, SO::ColumnMajor>, SWE::n_auxiliaries> aux_at_gp;

    DynMatrix<double> Fx_at_gp;
    DynMatrix<double> Fy_at_gp;

    DynMatrix<double> source_at_gp;
    DynMatrix<double> dbath_at_gp;
    DynMatrix<double> tau_s_at_gp;
    DynMatrix<double> dp_atm_at_gp;
    DynMatrix<double> dtide_pot_at_gp;

    HybMatrix<double, SWE::n_variables> q_prev_at_gp;
    HybMatrix<double, SWE::n_variables> del_q_DT_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> kronecker_DT_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dFx_dq_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dFy_dq_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dsource_dq_at_gp;

    DynMatrix<double> delta_local_inv;
    DynMatrix<double> delta_local;
    DynVector<double> rhs_local;
    DynVector<double> rhs_prev;
/*
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
*/
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
};
}

#endif
