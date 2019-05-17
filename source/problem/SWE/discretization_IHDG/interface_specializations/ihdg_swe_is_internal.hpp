#ifndef IHDG_SWE_IS_INTERNAL_HPP
#define IHDG_SWE_IS_INTERNAL_HPP

namespace SWE {
namespace IHDG {
namespace ISP {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface);

    template <typename EdgeInterfaceType>
    void ComputeInitTrace(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);
};

template <typename InterfaceType>
void Internal::Initialize(InterfaceType& intface) {} /*nothing to initialize*/

template <typename EdgeInterfaceType>
void Internal::ComputeInitTrace(EdgeInterfaceType& edge_int) {
    auto& intface = edge_int.interface;

    auto& state_in    = intface.data_in.state[0];
    auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

    auto& state_ex    = intface.data_ex.state[0];
    auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

    auto& edge_state    = edge_int.edge_data.edge_state;
    auto& edge_internal = edge_int.edge_data.edge_internal;

    for ( uint var = 0; var < SWE::n_variables; ++var ) {
        boundary_in.q_at_gp[var] = intface.ComputeUgpIN(state_in.q[var]);
        boundary_ex.q_at_gp[var] = intface.ComputeUgpEX(state_ex.q[var]);

        // Our definition of numerical flux implies q_hat = 0.5 * (q_in + q_ex)
        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            edge_internal.q_hat_at_gp(var, gp) =
                (boundary_in.q_at_gp[var][gp] + boundary_ex.q_at_gp[var][gp]) / 2.0;
        }
    }

    edge_state.q_hat = edge_int.L2Projection(edge_internal.q_hat_at_gp);
}

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(boundary_in.delta_global_kernel_at_gp, gp)    = column(boundary_in.dF_hat_dq_at_gp, gp);
        column(boundary_ex.delta_global_kernel_at_gp, gp_ex) = column(boundary_ex.dF_hat_dq_at_gp, gp_ex);

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = column(boundary_in.dF_hat_dq_hat_at_gp, gp);
        column(edge_internal.delta_hat_global_kernel_at_gp, gp) += column(boundary_ex.dF_hat_dq_hat_at_gp, gp_ex);

        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            edge_internal.rhs_global_kernel_at_gp(var, gp) = boundary_in.F_hat_at_gp[var][gp] + boundary_ex.F_hat_at_gp[var][gp];
        }
    }
}
}
}
}

#endif