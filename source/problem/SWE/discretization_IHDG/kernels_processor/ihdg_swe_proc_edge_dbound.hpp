#ifndef IHDG_SWE_PROC_EDGE_DBOUND_HPP
#define IHDG_SWE_PROC_EDGE_DBOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeDistributedType>
void Problem::local_edge_distributed_kernel(const RKStepper& stepper, EdgeDistributedType& edge_dbound) {
    auto& edge_state = edge_dbound.edge_data.edge_state;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    /* Take average of in/ex state as initial trace state */
    double ze_ex, qx_ex, qy_ex;
    double ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex;

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(
            gp, ze_ex, qx_ex, qy_ex, ze_flux_dot_n_ex, qx_flux_dot_n_ex, qy_flux_dot_n_ex);

        edge_state.ze_avg_at_gp[gp] = (boundary.ze_at_gp[gp] + ze_ex) / 2.0;
        edge_state.qx_avg_at_gp[gp] = (boundary.qx_at_gp[gp] + qx_ex) / 2.0;
        edge_state.qy_avg_at_gp[gp] = (boundary.qy_at_gp[gp] + qy_ex) / 2.0;
    }

    edge_dbound.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_dbound.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_dbound.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);
}
}
}

#endif