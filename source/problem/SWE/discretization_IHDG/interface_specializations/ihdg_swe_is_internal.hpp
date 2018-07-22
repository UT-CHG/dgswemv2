#ifndef IHDG_SWE_IS_INTERNAL_HPP
#define IHDG_SWE_IS_INTERNAL_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
namespace IS {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        boundary_in.delta_global_kernel_at_gp[gp]    = boundary_in.dF_hat_dq_at_gp[gp];
        boundary_ex.delta_global_kernel_at_gp[gp_ex] = boundary_ex.dF_hat_dq_at_gp[gp_ex];

        row(edge_internal.delta_hat_global_kernel_at_gp, gp) = boundary_in.dF_hat_dq_hat_at_gp[gp];
        row(edge_internal.delta_hat_global_kernel_at_gp, gp) += boundary_ex.dF_hat_dq_hat_at_gp[gp_ex];

        row(edge_internal.rhs_global_kernel_at_gp, gp) = row(boundary_in.F_hat_at_gp, gp);
        row(edge_internal.rhs_global_kernel_at_gp, gp) += row(boundary_ex.F_hat_at_gp, gp_ex);
    }
}
}
}
}

#endif