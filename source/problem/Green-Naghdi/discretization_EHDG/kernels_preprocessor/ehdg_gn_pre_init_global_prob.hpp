#ifndef EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP
#define EHDG_GN_PRE_INIT_GLOBAL_PROB_HPP

namespace GN {
namespace EHDG {
void Problem::initialize_global_dc_problem_serial(ProblemDiscretizationType& discretization,
                                                  uint& dc_global_dof_offset) {
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
        boundary_ex.w1_w1_hat.resize(GN::n_dimensions * edge_int.interface.data_ex.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());

        // Initialize w2 containers
        boundary_in.w2_w1_hat.resize(edge_int.interface.data_in.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());
        boundary_ex.w2_w1_hat.resize(edge_int.interface.data_ex.get_ndof(),
                                     GN::n_dimensions * edge_int.edge_data.get_ndof());

        // Initialize w1_hat containers
        edge_internal.w1_hat_w1_hat.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                           GN::n_dimensions * edge_int.edge_data.get_ndof());
        edge_internal.w1_hat_rhs.resize(GN::n_dimensions * edge_int.edge_data.get_ndof());

        boundary_in.w1_hat_w1.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     GN::n_dimensions * edge_int.interface.data_in.get_ndof());
        boundary_ex.w1_hat_w1.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     GN::n_dimensions * edge_int.interface.data_ex.get_ndof());

        boundary_in.w1_hat_w2.resize(GN::n_dimensions * edge_int.edge_data.get_ndof(),
                                     edge_int.interface.data_in.get_ndof());
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
        edge_internal.w1_hat_rhs.resize(GN::n_dimensions * edge_bound.edge_data.get_ndof());

        // Initialize w2 containers
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
                                                             uint& dc_global_dof_offset) {
    Problem::initialize_global_dc_problem_serial(discretization, dc_global_dof_offset);

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&dc_global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            // Set indexes for global matrix construction
            edge_internal.dc_global_dof_indx.resize(ndof_global * GN::n_dimensions);

            for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
                edge_internal.dc_global_dof_indx[indx] = indx + dc_global_dof_offset;
            }

            boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;

            dc_global_dof_offset += ndof_global * GN::n_dimensions;
        }

        // Initialize w1 containers
        boundary.w1_w1_hat.resize(GN::n_dimensions * edge_dbound.boundary.data.get_ndof(),
                                  GN::n_dimensions * edge_dbound.edge_data.get_ndof());

        // Initialize w2 containers
        boundary.w2_w1_hat.resize(edge_dbound.boundary.data.get_ndof(),
                                  GN::n_dimensions * edge_dbound.edge_data.get_ndof());

        // Initialize w1_hat containers
        edge_internal.w1_hat_w1_hat.resize(GN::n_dimensions * edge_dbound.edge_data.get_ndof(),
                                           GN::n_dimensions * edge_dbound.edge_data.get_ndof());
        edge_internal.w1_hat_rhs.resize(GN::n_dimensions * edge_dbound.edge_data.get_ndof());

        boundary.w1_hat_w1.resize(GN::n_dimensions * edge_dbound.edge_data.get_ndof(),
                                  GN::n_dimensions * edge_dbound.boundary.data.get_ndof());

        boundary.w1_hat_w2.resize(GN::n_dimensions * edge_dbound.edge_data.get_ndof(),
                                  edge_dbound.boundary.data.get_ndof());
    });
}

void Problem::initialize_global_dc_problem_parallel_finalize_pre_send(ProblemDiscretizationType& discretization,
                                                                      uint dc_global_dof_offset) {
    discretization.mesh_skeleton.CallForEachEdgeInterface([&dc_global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint ndof_global = edge_int.edge_data.get_ndof();

        // Set indexes for global matrix construction
        for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
            edge_internal.dc_global_dof_indx[indx] += dc_global_dof_offset;
        }

        boundary_in.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
        boundary_ex.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&dc_global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint ndof_global = edge_bound.edge_data.get_ndof();

        // Set indexes for global matrix construction
        for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
            edge_internal.dc_global_dof_indx[indx] += dc_global_dof_offset;
        }

        boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
    });

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&dc_global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            // Set indexes for global matrix construction
            for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
                edge_internal.dc_global_dof_indx[indx] += dc_global_dof_offset;
            }

            boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;

            std::vector<double> message(1);

            message[0] = (double)edge_internal.dc_global_dof_indx[0];

            edge_dbound.boundary.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dc_global_dof_indx, message);
        }
    });
}

template <typename Communicator>
void Problem::initialize_global_dc_problem_parallel_post_receive(ProblemDiscretizationType& discretization,
                                                                 Communicator& communicator,
                                                                 std::vector<uint>& dc_global_dof_indx) {
    discretization.mesh_skeleton.CallForEachEdgeInterface([&dc_global_dof_indx](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        edge_internal.dc_sol_offset = dc_global_dof_indx.size();

        dc_global_dof_indx.insert(
            dc_global_dof_indx.end(), edge_internal.dc_global_dof_indx.begin(), edge_internal.dc_global_dof_indx.end());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&dc_global_dof_indx](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        edge_internal.dc_sol_offset = dc_global_dof_indx.size();

        dc_global_dof_indx.insert(
            dc_global_dof_indx.end(), edge_internal.dc_global_dof_indx.begin(), edge_internal.dc_global_dof_indx.end());
    });

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&dc_global_dof_indx](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            std::vector<double> message(1);

            edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dc_global_dof_indx,
                                                                                   message);

            uint dc_global_dof_offset = (uint)message[0];

            edge_internal.dc_global_dof_indx.resize(ndof_global * GN::n_dimensions);

            // Set indexes for global matrix construction
            for (uint indx = 0; indx < ndof_global * GN::n_dimensions; ++indx) {
                edge_internal.dc_global_dof_indx[indx] = indx + dc_global_dof_offset;
            }

            boundary.dc_global_dof_indx = edge_internal.dc_global_dof_indx;
        }

        edge_internal.dc_sol_offset = dc_global_dof_indx.size();

        dc_global_dof_indx.insert(
            dc_global_dof_indx.end(), edge_internal.dc_global_dof_indx.begin(), edge_internal.dc_global_dof_indx.end());
    });
}
}
}

#endif