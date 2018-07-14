#ifndef IHDG_SWE_PROC_SOLVE_GLOBAL_PROB_HPP
#define IHDG_SWE_PROC_SOLVE_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename SimulationType>
bool Problem::solve_global_problem(SimulationType* simulation) {
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

    simulation->rhs_global =
        simulation->rhs_global - simulation->delta_global * simulation->delta_local_inv * simulation->rhs_local;

    int ipiv[simulation->rhs_global.size()];

    blaze::gesv(simulation->global, simulation->rhs_global, ipiv);

    simulation->rhs_local =
        simulation->delta_local_inv * (simulation->rhs_local - simulation->delta_hat_local * simulation->rhs_global);

    simulation->mesh.CallForEachElement([simulation](auto& elt) {
        const uint stage = simulation->stepper.GetStage();

        auto& state    = elt.data.state[stage + 1];
        auto& internal = elt.data.internal;

        uint elt_ID = elt.GetID();

        auto rhs_local = blaze::subvector(simulation->rhs_local, elt_ID * 9, 9);

        for (uint dof = 0; dof < elt.data.get_ndof(); ++dof) {
            state.q[dof][SWE::Variables::ze] += rhs_local[3 * dof];
            state.q[dof][SWE::Variables::qx] += rhs_local[3 * dof + 1];
            state.q[dof][SWE::Variables::qy] += rhs_local[3 * dof + 2];
        }
    });

    simulation->mesh_skeleton.CallForEachEdgeInterface([simulation](auto& edge_int) {
        auto& edge_state    = edge_int.edge_data.edge_state;
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint elt_in_ID = boundary_in.elt_ID;
        uint elt_ex_ID = boundary_ex.elt_ID;

        uint edg_ID = edge_int.GetID();

        auto rhs_global = blaze::subvector(simulation->rhs_global, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_int.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global[3 * dof + 2];
        }
    });

    simulation->mesh_skeleton.CallForEachEdgeBoundary([simulation](auto& edge_bound) {
        auto& edge_state    = edge_bound.edge_data.edge_state;
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint elt_ID = boundary.elt_ID;
        uint edg_ID = edge_bound.GetID();

        auto rhs_global = blaze::subvector(simulation->rhs_global, edg_ID * 6, 6);

        for (uint dof = 0; dof < edge_bound.edge_data.get_ndof(); ++dof) {
            edge_state.q_hat[dof][SWE::Variables::ze] += rhs_global[3 * dof];
            edge_state.q_hat[dof][SWE::Variables::qx] += rhs_global[3 * dof + 1];
            edge_state.q_hat[dof][SWE::Variables::qy] += rhs_global[3 * dof + 2];
        }
    });

    double del_norm = std::hypot(blaze::norm(simulation->rhs_local_prev - simulation->rhs_local),
                                 blaze::norm(simulation->rhs_global_prev - simulation->rhs_global));

    double norm = std::hypot(blaze::norm(simulation->rhs_local), blaze::norm(simulation->rhs_global));

    simulation->rhs_local_prev  = simulation->rhs_local;
    simulation->rhs_global_prev = simulation->rhs_global;

    if (del_norm / norm < 1e-3) {
        return true;
    }

    return false;
}
}
}

#endif