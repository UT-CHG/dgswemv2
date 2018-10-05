#ifndef EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP
#define EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_global_problem(ProblemDiscretizationType& discretization) {
    SWE::EHDG::Problem::initialize_global_problem(discretization);

    uint local_dof_offset  = 0;
    uint global_dof_offset = 0;

    discretization.mesh.CallForEachElement([&local_dof_offset](auto& elt) {
        auto& internal = elt.data.internal;

        // Set offsets for global matrix construction
        internal.local_dof_offset = local_dof_offset;
        local_dof_offset += elt.data.get_ndof();

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].local_dof_offset = internal.local_dof_offset;
        }

        // Initialize w1 containers
        internal.w1_w1.resize(GN::n_dimensions * elt.data.get_ndof(), GN::n_dimensions * elt.data.get_ndof());
        internal.w1_w2.resize(GN::n_dimensions * elt.data.get_ndof(), elt.data.get_ndof());
        internal.w1_rhs.resize(GN::n_dimensions * elt.data.get_ndof());

        // Initialize w2 containers
        internal.w2_w1.resize(elt.data.get_ndof(), GN::n_dimensions * elt.data.get_ndof());
        internal.w2_w2.resize(elt.data.get_ndof(), elt.data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;

        boundary_in.global_dof_offset = global_dof_offset;
        boundary_ex.global_dof_offset = global_dof_offset;

        global_dof_offset += edge_int.edge_data.get_ndof();

        // Initialize w1 containers
        boundary_in.w1_w1_hat.resize(GN::n_dimensions * edge_int.interface.data_in.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());
        boundary_in.w2_w1_hat.resize(edge_int.interface.data_in.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());

        // Initialize w2 containers
        boundary_ex.w1_w1_hat.resize(GN::n_dimensions * edge_int.interface.data_ex.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());
        boundary_ex.w2_w1_hat.resize(edge_int.interface.data_ex.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());

        // Initialize w1_hat containers
        edge_internal.w1_hat_w1_hat.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                           GN::n_dimensions * edge_int.edge_data.get_ndof());

        boundary_in.w1_hat_w1.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     GN::n_dimensions * edge_int.interface.data_in.get_ndof());

        boundary_in.w1_hat_w2.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     edge_int.interface.data_in.get_ndof());

        boundary_ex.w1_hat_w1.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     GN::n_dimensions * edge_int.interface.data_ex.get_ndof());

        boundary_ex.w1_hat_w2.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     edge_int.interface.data_ex.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;

        boundary.global_dof_offset = global_dof_offset;

        global_dof_offset += edge_bound.edge_data.get_ndof();

        // Initialize w1 containers
        boundary.w1_w1_hat.resize(GN::n_dimensions * edge_bound.boundary.data.get_ndof(),
                                  GN::n_dimensions * edge_bound.edge_data.get_ndof());
        boundary.w2_w1_hat.resize(edge_bound.boundary.data.get_ndof(),
                                  GN::n_dimensions * edge_bound.edge_data.get_ndof());

        // Initialize w1_hat containers
        edge_internal.w1_hat_w1_hat.resize(GN::n_dimensions * edge_bound.edge_data.get_ndof(),
                                           GN::n_dimensions * edge_bound.edge_data.get_ndof());

        boundary.w1_hat_w1.resize(GN::n_dimensions * edge_bound.edge_data.get_ndof(),
                                  GN::n_dimensions * edge_bound.boundary.data.get_ndof());
    });

    discretization.global_data.w1_w1_hat.resize(local_dof_offset * GN::n_dimensions,
                                                global_dof_offset * GN::n_dimensions);
    discretization.global_data.w1_rhs.resize(local_dof_offset * GN::n_dimensions);

    discretization.global_data.w2_w1.resize(local_dof_offset, local_dof_offset * GN::n_dimensions);
    discretization.global_data.w2_w2_inv.resize(local_dof_offset, local_dof_offset);
    discretization.global_data.w2_w1_hat.resize(local_dof_offset, global_dof_offset * GN::n_dimensions);

    discretization.global_data.w1_hat_w1.resize(global_dof_offset * GN::n_dimensions,
                                                local_dof_offset * GN::n_dimensions);
    discretization.global_data.w1_hat_w2.resize(global_dof_offset * GN::n_dimensions, local_dof_offset);
    discretization.global_data.w1_hat_w1_hat.resize(global_dof_offset * GN::n_dimensions,
                                                    global_dof_offset * GN::n_dimensions);
    discretization.global_data.w1_hat_rhs.resize(global_dof_offset * GN::n_dimensions);
}
}
}

#endif