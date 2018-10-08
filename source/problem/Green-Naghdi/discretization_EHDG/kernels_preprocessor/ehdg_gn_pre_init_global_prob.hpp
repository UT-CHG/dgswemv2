#ifndef EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP
#define EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_global_dc_problem(ProblemDiscretizationType& discretization, uint& dc_global_dof_offset) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        // Initialize w1 containers
        internal.w1_w1.resize(GN::n_dimensions * elt.data.get_ndof(), GN::n_dimensions * elt.data.get_ndof());
        internal.w1_w2.resize(GN::n_dimensions * elt.data.get_ndof(), elt.data.get_ndof());
        internal.w1_rhs.resize(GN::n_dimensions * elt.data.get_ndof());

        // Initialize w2 containers
        internal.w2_w1.resize(elt.data.get_ndof(), GN::n_dimensions * elt.data.get_ndof());
        internal.w2_w2.resize(elt.data.get_ndof(), elt.data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&dc_global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint ndof_global = edge_int.edge_data.get_ndof();

        // Set indexes for global matrix construction
        edge_internal.dc_global_dof_indx.resize(ndof_global * GN::n_dimensions);

        for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
            edge_internal.dc_global_dof_indx[indx] = indx + dc_global_dof_offset;
        }

        boundary_in.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
        boundary_ex.dc_global_dof_indx = edge_internal.dc_global_dof_indx;

        dc_global_dof_offset += ndof_global * GN::n_dimensions;

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

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&dc_global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint ndof_global = edge_bound.edge_data.get_ndof();

        // Set indexes for global matrix construction
        edge_internal.dc_global_dof_indx.resize(ndof_global * GN::n_dimensions);

        for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
            edge_internal.dc_global_dof_indx[indx] = indx + dc_global_dof_offset;
        }

        boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;

        dc_global_dof_offset += ndof_global * GN::n_dimensions;

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
}

template <typename Communicator>
void Problem::initialize_global_dc_problem_parallel_pre_send(ProblemDiscretizationType& discretization,
                                                             Communicator& communicator,
                                                             uint& dc_global_dof_offset) {}

void Problem::initialize_global_dc_problem_parallel_finalize_pre_send(ProblemDiscretizationType& discretization,
                                                                      uint dc_global_dof_offset) {}

template <typename Communicator>
void Problem::initialize_global_dc_problem_parallel_post_receive(ProblemDiscretizationType& discretization,
                                                                 Communicator& communicator,
                                                                 std::vector<uint>& dc_global_dof_indx) {}
}
}

#endif