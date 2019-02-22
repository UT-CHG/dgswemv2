#ifndef SWE_DATA_BOUNDARY_HPP
#define SWE_DATA_BOUNDARY_HPP

namespace SWE {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : q_at_gp(SWE::n_variables, ngp),
          aux_at_gp(SWE::n_auxiliaries, ngp),
          F_hat_at_gp(SWE::n_variables, ngp),
          dF_hat_dq_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          dF_hat_dq_hat_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          delta_global_kernel_at_gp(SWE::n_variables * SWE::n_variables, ngp) {}

    HybMatrix<double, SWE::n_variables> q_at_gp;
    HybMatrix<double, SWE::n_auxiliaries> aux_at_gp;

    HybMatrix<double, SWE::n_variables> F_hat_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dF_hat_dq_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> dF_hat_dq_hat_at_gp;
    HybMatrix<double, SWE::n_variables * SWE::n_variables> delta_global_kernel_at_gp;

    DynMatrix<double> delta_global;
    DynMatrix<double> delta_hat_local;

    std::vector<uint> global_dof_indx;

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

#endif
