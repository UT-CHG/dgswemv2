#ifndef IHDG_GN_IS_INTERNAL_HPP
#define IHDG_GN_IS_INTERNAL_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/Green-Naghdi/discretization_IHDG/stabilization_parameters/IHDG_gn_stabilization_params.hpp"

namespace GN {
namespace IHDG {
namespace IS {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernelsSWE(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFluxSWE(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernelsDC(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernelsSWE(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(edge_internal.rhs_global_kernel_at_gp, gp) = column(boundary_in.Fn_at_gp, gp);
        column(edge_internal.rhs_global_kernel_at_gp, gp) += column(boundary_ex.Fn_at_gp, gp_ex);
    }

    // Add tau terms
    add_kernel_tau_terms_intface_LF(edge_int);
}

template <typename EdgeInterfaceType>
void Internal::ComputeNumericalFluxSWE(EdgeInterfaceType& edge_int) {
    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    boundary_in.F_hat_at_gp = boundary_in.Fn_at_gp;
    boundary_ex.F_hat_at_gp = boundary_ex.Fn_at_gp;

    // Add tau terms
    add_F_hat_tau_terms_intface_LF(edge_int);
}

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernelsDC(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double tau = -20;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(edge_internal.w1_hat_w1_hat_kernel_at_gp, gp) = -2.0 * tau * IdentityVector<double>(GN::n_dimensions);

        column(boundary_in.w1_hat_w1_kernel_at_gp, gp)    = tau * IdentityVector<double>(GN::n_dimensions);
        column(boundary_ex.w1_hat_w1_kernel_at_gp, gp_ex) = tau * IdentityVector<double>(GN::n_dimensions);
    }

    boundary_in.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_in;
    boundary_ex.w1_hat_w2_kernel_at_gp = edge_int.interface.surface_normal_ex;
}
}
}
}

#endif