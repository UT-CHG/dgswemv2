#ifndef EHDG_SWE_IS_INTERNAL_HPP
#define EHDG_SWE_IS_INTERNAL_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_EHDG/stabilization_parameters/ehdg_swe_stabilization_params.hpp"

namespace SWE {
namespace EHDG {
namespace IS {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFlux(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    auto& edge_state  = edge_int.edge_data.edge_state;
    auto& edge_global = edge_int.edge_data.edge_global;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        // Add F_in * n_in terms
        edge_global.ze_rhs_kernel_at_gp[gp] = -boundary_in.ze_flux_dot_n_at_gp[gp];
        edge_global.qx_rhs_kernel_at_gp[gp] = -boundary_in.qx_flux_dot_n_at_gp[gp];
        edge_global.qy_rhs_kernel_at_gp[gp] = -boundary_in.qy_flux_dot_n_at_gp[gp];

        // Add F_ex * n_ex terms
        edge_global.ze_rhs_kernel_at_gp[gp] += -boundary_ex.ze_flux_dot_n_at_gp[gp_ex];
        edge_global.qx_rhs_kernel_at_gp[gp] += -boundary_ex.qx_flux_dot_n_at_gp[gp_ex];
        edge_global.qy_rhs_kernel_at_gp[gp] += -boundary_ex.qy_flux_dot_n_at_gp[gp_ex];
    }

    // Add tau terms
    add_kernel_tau_terms_intface_LF(edge_int);
}

template <typename EdgeInterfaceType>
void Internal::ComputeNumericalFlux(EdgeInterfaceType& edge_int) {
    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        boundary_in.ze_numerical_flux_at_gp[gp] = boundary_in.ze_flux_dot_n_at_gp[gp];
        boundary_in.qx_numerical_flux_at_gp[gp] = boundary_in.qx_flux_dot_n_at_gp[gp];
        boundary_in.qy_numerical_flux_at_gp[gp] = boundary_in.qy_flux_dot_n_at_gp[gp];

        boundary_ex.ze_numerical_flux_at_gp[gp_ex] = boundary_ex.ze_flux_dot_n_at_gp[gp_ex];
        boundary_ex.qx_numerical_flux_at_gp[gp_ex] = boundary_ex.qx_flux_dot_n_at_gp[gp_ex];
        boundary_ex.qy_numerical_flux_at_gp[gp_ex] = boundary_ex.qy_flux_dot_n_at_gp[gp_ex];
    }

    // Add tau terms
    add_flux_tau_terms_intface_LF(edge_int);
}
}
}
}

#endif