#ifndef IHDG_SWE_PROC_ASM_GLOBAL_PROB_HPP
#define IHDG_SWE_PROC_ASM_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename SimulationType>
void Problem::assemble_global_problem(SimulationType* simulation) {
    simulation->mesh.CallForEachElement([simulation](auto& elt) {
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        blaze::submatrix(simulation->delta_local_inv, elt_ID * 9, elt_ID * 9, 9, 9) = inv(internal.delta_local);
        blaze::subvector(simulation->rhs_local, elt_ID * 9, 9)                      = internal.rhs_local;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); bound_id++) {
            uint edg_ID = elt.data.boundary[bound_id].edg_ID;

            blaze::submatrix(simulation->delta_hat_local, elt_ID * 9, edg_ID * 6, 9, 6) =
                elt.data.boundary[bound_id].delta_hat_local;
        }
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface([simulation](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        blaze::submatrix(simulation->delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        blaze::subvector(simulation->rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        blaze::submatrix(simulation->delta_global, edg_ID * 6, elt_in_ID * 9, 6, 9) = boundary_in.delta_global;
        blaze::submatrix(simulation->delta_global, edg_ID * 6, elt_ex_ID * 9, 6, 9) = boundary_ex.delta_global;
    });

    simulation->mesh_skeleton.CallForEachEdgeBoundary([simulation](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        blaze::submatrix(simulation->delta_hat_global, edg_ID * 6, edg_ID * 6, 6, 6) = edge_internal.delta_hat_global;
        blaze::subvector(simulation->rhs_global, edg_ID * 6, 6)                      = edge_internal.rhs_global;

        blaze::submatrix(simulation->delta_global, edg_ID * 6, elt_ID * 9, 6, 9) = boundary.delta_global;
    });

    simulation->global = simulation->delta_hat_global -
                         simulation->delta_global * simulation->delta_local_inv * simulation->delta_hat_local;

    simulation->rhs =
        simulation->rhs_global - simulation->delta_global * simulation->delta_local_inv * simulation->rhs_local;
}
}
}

#endif