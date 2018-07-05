#ifndef IHDG_SWE_PROC_EDGE_BOUND_HPP
#define IHDG_SWE_PROC_EDGE_BOUND_HPP

#include <eigen3/Eigen/Dense>
#include "problem/SWE/discretization_IHDG/stabilization_parameters/ihdg_swe_stabilization_params.hpp"

namespace SWE {
namespace IHDG {
template <typename EdgeBoundaryType>
void Problem::prepare_edge_boundary_kernel(const RKStepper& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_state = edge_bound.edge_data.edge_state;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    /* Take in state as initial edge state */
    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        edge_state.ze_avg_at_gp[gp] = boundary.ze_at_gp[gp];
        edge_state.qx_avg_at_gp[gp] = boundary.qx_at_gp[gp];
        edge_state.qy_avg_at_gp[gp] = boundary.qy_at_gp[gp];
    }

    edge_bound.L2Projection(edge_state.ze_avg_at_gp, edge_state.ze_hat);
    edge_bound.L2Projection(edge_state.qx_avg_at_gp, edge_state.qx_hat);
    edge_bound.L2Projection(edge_state.qy_avg_at_gp, edge_state.qy_hat);
}
}
}

#endif