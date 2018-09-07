#ifndef IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
void Problem::initialize_global_problem_serial(HDGDiscretization<Problem>& discretization) {
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
void Problem::initialize_global_problem_parallel_pre_send(HDGDiscretization<Problem>& discretization,
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

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        set_constant(edge_dbound.edge_data.edge_state.q_hat, 0.0);

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            edge_internal.global_dof_offset = global_dof_offset;

            boundary.global_dof_offset = global_dof_offset;

            global_dof_offset += edge_dbound.edge_data.get_ndof();
        } else {
            edge_internal.global_dof_offset = 0;

            boundary.global_dof_offset = 0;
        }

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

void Problem::initialize_global_problem_parallel_finalize_pre_send(HDGDiscretization<Problem>& discretization,
                                                                   uint local_dof_offset,
                                                                   uint global_dof_offset) {
    discretization.mesh.CallForEachElement([&local_dof_offset](auto& elt) {
        auto& internal = elt.data.internal;

        // Set offsets for global matrix construction
        internal.local_dof_offset += local_dof_offset;

        for (uint bound_id = 0; bound_id < elt.data.get_nbound(); ++bound_id) {
            elt.data.boundary[bound_id].local_dof_offset += local_dof_offset;
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
void Problem::initialize_global_problem_parallel_post_receive(HDGDiscretization<Problem>& discretization,
                                                              Communicator& communicator) {
    discretization.mesh_skeleton.CallForEachEdgeDistributed([](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        std::vector<double> message(1);

        edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(SWE::CommTypes::preprocessor, message);

        uint global_dof_offset = (uint)message[0];

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            edge_internal.global_dof_offset = global_dof_offset;

            boundary.global_dof_offset = global_dof_offset;
        }
    });
}
}
}

#endif