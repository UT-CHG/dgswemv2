#ifndef IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
void Problem::initialize_global_problem_serial_kernel(HDGDiscretization<Problem>& discretization) {
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

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        set_constant(edge_int.edge_data.edge_state.q_hat, 0.0);

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;

        boundary_in.global_dof_offset = global_dof_offset;
        boundary_ex.global_dof_offset = global_dof_offset;

        global_dof_offset += edge_int.edge_data.get_ndof();

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                              SWE::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_hat_local containers
        boundary_in.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_in.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());
        boundary_ex.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_ex.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_global containers
        boundary_in.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_in.get_ndof());
        boundary_ex.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_ex.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        set_constant(edge_bound.edge_data.edge_state.q_hat, 0.0);

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;

        boundary.global_dof_offset = global_dof_offset;

        global_dof_offset += edge_bound.edge_data.get_ndof();

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_hat_local container
        boundary.delta_hat_local.resize(SWE::n_variables * edge_bound.boundary.data.get_ndof(),
                                        SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_global container
        boundary.delta_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                     SWE::n_variables * edge_bound.boundary.data.get_ndof());
    });

    discretization.global_data.delta_local_inv.resize(local_dof_offset * SWE::n_variables,
                                                      local_dof_offset * SWE::n_variables);
    discretization.global_data.delta_hat_local.resize(local_dof_offset * SWE::n_variables,
                                                      global_dof_offset * SWE::n_variables);
    discretization.global_data.rhs_local.resize(local_dof_offset * SWE::n_variables);

    discretization.global_data.delta_global.resize(global_dof_offset * SWE::n_variables,
                                                   local_dof_offset * SWE::n_variables);
    discretization.global_data.delta_hat_global.resize(global_dof_offset * SWE::n_variables,
                                                       global_dof_offset * SWE::n_variables);
    discretization.global_data.rhs_global.resize(global_dof_offset * SWE::n_variables);
}

template <typename Communicator>
void Problem::initialize_global_problem_parallel_pre_send_kernel(HDGDiscretization<Problem>& discretization,
                                                                 Communicator& communicator,
                                                                 uint& local_dof_offset,
                                                                 uint& global_dof_offset) {
    discretization.mesh.CallForEachElement([&local_dof_offset](auto& elt) {
        auto& internal = elt.data.internal;

        // Set offsets for global matrix construction
        internal.local_dof_offset = local_dof_offset;
        local_dof_offset += elt.data.get_ndof();

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].local_dof_offset = internal.local_dof_offset;
        }

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        set_constant(edge_int.edge_data.edge_state.q_hat, 0.0);

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;

        boundary_in.global_dof_offset = global_dof_offset;
        boundary_ex.global_dof_offset = global_dof_offset;

        global_dof_offset += edge_int.edge_data.get_ndof();

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                              SWE::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_hat_local containers
        boundary_in.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_in.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());
        boundary_ex.delta_hat_local.resize(SWE::n_variables * edge_int.interface.data_ex.get_ndof(),
                                           SWE::n_variables * edge_int.edge_data.get_ndof());

        // Initialize delta_global containers
        boundary_in.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_in.get_ndof());
        boundary_ex.delta_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                        SWE::n_variables * edge_int.interface.data_ex.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        set_constant(edge_bound.edge_data.edge_state.q_hat, 0.0);

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset = global_dof_offset;

        boundary.global_dof_offset = global_dof_offset;

        global_dof_offset += edge_bound.edge_data.get_ndof();

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_hat_local container
        boundary.delta_hat_local.resize(SWE::n_variables * edge_bound.boundary.data.get_ndof(),
                                        SWE::n_variables * edge_bound.edge_data.get_ndof());

        // Initialize delta_global container
        boundary.delta_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                     SWE::n_variables * edge_bound.boundary.data.get_ndof());
    });

    std::vector<bool> set_global_dof_offset;

    // set who numbers the global dofs (lowest locality/submesh side)
    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); ++rank_boundary_id) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        uint n_dbound = rank_boundary.db_data.elements_in.size();

        if (rank_boundary.db_data.locality_in < rank_boundary.db_data.locality_ex ||
            (rank_boundary.db_data.locality_in == rank_boundary.db_data.locality_ex &&
             rank_boundary.db_data.submesh_in < rank_boundary.db_data.submesh_ex)) {
            // *** //
            set_global_dof_offset.insert(set_global_dof_offset.begin(), n_dbound, true);
        } else {
            set_global_dof_offset.insert(set_global_dof_offset.begin(), n_dbound, false);
        }
    }

    discretization.mesh_skeleton.CallForEachEdgeDistributed(
        [&global_dof_offset, &set_global_dof_offset](auto& edge_dbound) {
            auto& edge_internal = edge_dbound.edge_data.edge_internal;

            auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

            set_constant(edge_dbound.edge_data.edge_state.q_hat, 0.0);

            if (set_global_dof_offset.back()) {
                edge_internal.global_dof_offset = global_dof_offset;

                boundary.global_dof_offset = global_dof_offset;

                global_dof_offset += edge_dbound.edge_data.get_ndof();
            } else {
                edge_internal.global_dof_offset = 0;

                boundary.global_dof_offset = 0;
            }

            set_global_dof_offset.pop_back();

            // Initialize delta_hat_global and rhs_global containers
            edge_internal.delta_hat_global.resize(SWE::n_variables * edge_dbound.edge_data.get_ndof(),
                                                  SWE::n_variables * edge_dbound.edge_data.get_ndof());
            edge_internal.rhs_global.resize(SWE::n_variables * edge_dbound.edge_data.get_ndof());

            // Initialize delta_hat_local container
            boundary.delta_hat_local.resize(SWE::n_variables * edge_dbound.boundary.data.get_ndof(),
                                            SWE::n_variables * edge_dbound.edge_data.get_ndof());

            // Initialize delta_global container
            boundary.delta_global.resize(SWE::n_variables * edge_dbound.edge_data.get_ndof(),
                                         SWE::n_variables * edge_dbound.boundary.data.get_ndof());
        });
}

void Problem::initialize_global_problem_parallel_finalize_pre_send_kernel(HDGDiscretization<Problem>& discretization,
                                                                          uint local_dof_offset,
                                                                          uint global_dof_offset) {
    discretization.mesh.CallForEachElement([&local_dof_offset](auto& elt) {
        auto& internal = elt.data.internal;

        // Set offsets for global matrix construction
        internal.local_dof_offset += local_dof_offset;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].local_dof_offset += internal.local_dof_offset;
        }
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset += global_dof_offset;

        boundary_in.global_dof_offset += global_dof_offset;
        boundary_ex.global_dof_offset += global_dof_offset;
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset += global_dof_offset;

        boundary.global_dof_offset += global_dof_offset;
    });

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        // Set offsets for global matrix construction
        edge_internal.global_dof_offset += global_dof_offset;

        boundary.global_dof_offset += global_dof_offset;

        std::vector<double> message(1);

        message[0] = (double)edge_internal.global_dof_offset;

        edge_dbound.boundary.boundary_condition.exchanger.SetToSendBuffer(SWE::CommTypes::preprocessor, message);
    });
}

template <typename Communicator>
void Problem::initialize_global_problem_parallel_post_receive_kernel(HDGDiscretization<Problem>& discretization,
                                                                     Communicator& communicator) {
    std::vector<bool> set_global_dof_offset;

    // set who numbers the global dofs (lowest locality/submesh side)
    for (uint rank_boundary_id = 0; rank_boundary_id < communicator.GetRankBoundaryNumber(); ++rank_boundary_id) {
        typename Communicator::RankBoundaryType& rank_boundary = communicator.GetRankBoundary(rank_boundary_id);

        uint n_dbound = rank_boundary.db_data.elements_in.size();

        if (rank_boundary.db_data.locality_in < rank_boundary.db_data.locality_ex ||
            (rank_boundary.db_data.locality_in == rank_boundary.db_data.locality_ex &&
             rank_boundary.db_data.submesh_in < rank_boundary.db_data.submesh_ex)) {
            // *** //
            set_global_dof_offset.insert(set_global_dof_offset.begin(), n_dbound, true);
        } else {
            set_global_dof_offset.insert(set_global_dof_offset.begin(), n_dbound, false);
        }
    }

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&set_global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        std::vector<double> message(1);

        edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(SWE::CommTypes::preprocessor, message);

        uint global_dof_offset = (uint)message[0];

        if (!set_global_dof_offset.back()) {
            edge_internal.global_dof_offset = global_dof_offset;

            boundary.global_dof_offset = global_dof_offset;
        }

        set_global_dof_offset.pop_back();
    });
}
}
}

#endif