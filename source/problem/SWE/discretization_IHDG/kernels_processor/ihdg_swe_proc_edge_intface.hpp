#ifndef IHDG_SWE_PROC_EDGE_INTFACE_HPP
#define IHDG_SWE_PROC_EDGE_INTFACE_HPP

#include <eigen3/Eigen/Dense>

namespace SWE {
namespace IHDG {
template <typename EdgeInterfaceType>
void Problem::prepare_edge_interface_kernel(const RKStepper& stepper, EdgeInterfaceType& edge_int) {
    auto& edge_state = edge_int.edge_data.edge_state;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    /* Take average of in/ex state as initial trace state */
    uint gp_ex = 0;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        edge_state.ze_avg_at_gp[gp] = (boundary_in.ze_at_gp[gp] + boundary_ex.ze_at_gp[gp_ex]) / 2.0;
        edge_state.qx_avg_at_gp[gp] = (boundary_in.qx_at_gp[gp] + boundary_ex.qx_at_gp[gp_ex]) / 2.0;
        edge_state.qy_avg_at_gp[gp] = (boundary_in.qy_at_gp[gp] + boundary_ex.qy_at_gp[gp_ex]) / 2.0;
    }

    edge_int.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_int.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_int.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);
}
}
}

#endif