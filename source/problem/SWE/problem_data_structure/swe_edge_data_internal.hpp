#ifndef SWE_EDGE_DATA_INTERNAL_HPP
#define SWE_EDGE_DATA_INTERNAL_HPP

namespace SWE {
struct EdgeInternal {
    EdgeInternal() = default;
    EdgeInternal(const uint ngp)
        : q_hat_at_gp(SWE::n_variables, ngp),
          aux_hat_at_gp(SWE::n_auxiliaries, ngp),
          q_init_at_gp(SWE::n_variables, ngp),
          tau(ngp),
          dtau_dze(ngp),
          dtau_dqx(ngp),
          dtau_dqy(ngp),
          delta_hat_global_kernel_at_gp(SWE::n_variables * SWE::n_variables, ngp),
          rhs_global_kernel_at_gp(SWE::n_variables, ngp) {}

    HybMatrix<double, SWE::n_variables> q_hat_at_gp;
    HybMatrix<double, SWE::n_auxiliaries> aux_hat_at_gp;
    HybMatrix<double, SWE::n_variables> q_init_at_gp;

    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> tau;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dtau_dze;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dtau_dqx;
    AlignedVector<StatMatrix<double, SWE::n_variables, SWE::n_variables>> dtau_dqy;

    HybMatrix<double, SWE::n_variables * SWE::n_variables> delta_hat_global_kernel_at_gp;
    HybMatrix<double, SWE::n_variables> rhs_global_kernel_at_gp;

    DynMatrix<double> delta_hat_global;
    DynVector<double> rhs_global;

    DynVector<double> delta_hat_global_flat;
    std::vector<DynVector<double>> delta_hat_global_con_flat;

    std::vector<uint> global_dof_indx;
    uint sol_offset;
};
}

#endif
