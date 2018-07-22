#ifndef IHDG_SWE_PROC_ASM_GLOBAL_PROB_HPP
#define IHDG_SWE_PROC_ASM_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename SimulationType>
void Problem::initialize_iteration(SimulationType* simulation) {
    simulation->mesh.CallForEachElement([simulation](auto& elt) {
        const uint stage = simulation->stepper.GetStage();

        auto& state      = elt.data.state[stage + 1];
        auto& state_prev = elt.data.state[stage];
        auto& internal   = elt.data.internal;

        state.q = state_prev.q;

        internal.q_prev_at_gp = elt.ComputeUgp(state_prev.q);
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface([simulation](auto& edge_int) {
        const uint stage = simulation->stepper.GetStage();

        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& state_in    = edge_int.interface.data_in.state[stage + 1];
        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];

        auto& state_ex    = edge_int.interface.data_ex.state[stage + 1];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        boundary_in.q_at_gp = edge_int.interface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = edge_int.interface.ComputeUgpEX(state_ex.q);

        /* Take average of in/ex state as initial trace state */
        uint gp_ex;
        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

            column(edge_internal.q_init_at_gp, gp) =
                (column(boundary_in.q_at_gp, gp) + column(boundary_ex.q_at_gp, gp_ex)) / 2.0;
        }

        edge_state.q_hat = edge_int.L2Projection(edge_internal.q_init_at_gp);
    });

    simulation->mesh_skeleton.CallForEachEdgeBoundary([simulation](auto& edge_bound) {
        const uint stage = simulation->stepper.GetStage();

        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& state    = edge_bound.boundary.data.state[stage + 1];
        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        boundary.q_at_gp = edge_bound.boundary.ComputeUgp(state.q);

        /* Take in state as initial edge state */
        edge_internal.q_init_at_gp = boundary.q_at_gp;

        edge_state.q_hat = edge_bound.L2Projection(edge_internal.q_init_at_gp);
    });
}
}
}

#endif