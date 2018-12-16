#ifndef EHDG_SWE_IS_INTERNAL_HPP
#define EHDG_SWE_IS_INTERNAL_HPP

namespace SWE {
namespace EHDG {
namespace ISP {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFlux(EdgeInterfaceType& edge_int);
};

template <typename InterfaceType>
void Internal::Initialize(InterfaceType& intface) {}

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    SWE::get_tau_LF(edge_internal.q_hat_at_gp,
                    edge_internal.aux_hat_at_gp,
                    edge_int.interface.surface_normal_in,
                    edge_internal.tau);
    SWE::get_dtau_dze_LF(edge_internal.q_hat_at_gp,
                         edge_internal.aux_hat_at_gp,
                         edge_int.interface.surface_normal_in,
                         edge_internal.dtau_dze);
    SWE::get_dtau_dqx_LF(edge_internal.q_hat_at_gp,
                         edge_internal.aux_hat_at_gp,
                         edge_int.interface.surface_normal_in,
                         edge_internal.dtau_dqx);
    SWE::get_dtau_dqy_LF(edge_internal.q_hat_at_gp,
                         edge_internal.aux_hat_at_gp,
                         edge_int.interface.surface_normal_in,
                         edge_internal.dtau_dqy);

    StatVector<double, SWE::n_variables> del_q;
    StatMatrix<double, SWE::n_variables, SWE::n_variables> dtau_delq;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        del_q = column(boundary_in.q_at_gp, gp) + column(boundary_ex.q_at_gp, gp_ex) -
                2.0 * column(edge_internal.q_hat_at_gp, gp);

        column(dtau_delq, SWE::Variables::ze) = edge_internal.dtau_dze[gp] * del_q;
        column(dtau_delq, SWE::Variables::qx) = edge_internal.dtau_dqx[gp] * del_q;
        column(dtau_delq, SWE::Variables::qy) = edge_internal.dtau_dqy[gp] * del_q;

        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            column(boundary_in.Fn_at_gp, gp) + column(boundary_ex.Fn_at_gp, gp_ex) + edge_internal.tau[gp] * del_q;

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) =
            flatten<double>(dtau_delq - 2.0 * edge_internal.tau[gp]);
    }
}

template <typename EdgeInterfaceType>
void Internal::ComputeNumericalFlux(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    SWE::get_tau_LF(edge_internal.q_hat_at_gp,
                    edge_internal.aux_hat_at_gp,
                    edge_int.interface.surface_normal_in,
                    edge_internal.tau);

    boundary_in.F_hat_at_gp = boundary_in.Fn_at_gp;
    boundary_ex.F_hat_at_gp = boundary_ex.Fn_at_gp;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(boundary_in.F_hat_at_gp, gp) +=
            edge_internal.tau[gp] * (column(boundary_in.q_at_gp, gp) - column(edge_internal.q_hat_at_gp, gp));
        column(boundary_ex.F_hat_at_gp, gp_ex) +=
            edge_internal.tau[gp] * (column(boundary_ex.q_at_gp, gp_ex) - column(edge_internal.q_hat_at_gp, gp));
    }
}
}
}
}

#endif