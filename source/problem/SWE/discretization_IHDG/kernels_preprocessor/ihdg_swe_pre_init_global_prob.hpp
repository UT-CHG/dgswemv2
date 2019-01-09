#ifndef IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename ProblemType>
void Problem::initialize_global_problem_serial(HDGDiscretization<ProblemType>& discretization,
                                               uint& global_dof_offset) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& internal = elt.data.internal;

        // Initialize delta_local and rhs_local containers
        internal.delta_local.resize(SWE::n_variables * elt.data.get_ndof(), SWE::n_variables * elt.data.get_ndof());
        internal.rhs_local.resize(SWE::n_variables * elt.data.get_ndof());
        internal.rhs_prev.resize(SWE::n_variables * elt.data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint ndof_global = edge_int.edge_data.get_ndof();

        // Set indexes for global matrix construction
        edge_internal.global_dof_indx.resize(ndof_global * SWE::n_variables);

        for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
            edge_internal.global_dof_indx[indx] = indx + global_dof_offset;
        }

        boundary_in.global_dof_indx = edge_internal.global_dof_indx;
        boundary_ex.global_dof_indx = edge_internal.global_dof_indx;

        global_dof_offset += ndof_global * SWE::n_variables;

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

        edge_int.interface.specialization.ComputeInitTrace(edge_int);
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint ndof_global = edge_bound.edge_data.get_ndof();

        // Set indexes for global matrix construction
        edge_internal.global_dof_indx.resize(ndof_global * SWE::n_variables);

        for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
            edge_internal.global_dof_indx[indx] = indx + global_dof_offset;
        }

        boundary.global_dof_indx = edge_internal.global_dof_indx;

        global_dof_offset += ndof_global * SWE::n_variables;

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

        edge_bound.boundary.boundary_condition.ComputeInitTrace(edge_bound);
    });
}

template <typename ProblemType, typename Communicator>
void Problem::initialize_global_problem_parallel_pre_send(HDGDiscretization<ProblemType>& discretization,
                                                          Communicator& communicator,
                                                          uint& global_dof_offset) {
    Problem::initialize_global_problem_serial(discretization, global_dof_offset);

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            // Set indexes for global matrix construction
            edge_internal.global_dof_indx.resize(ndof_global * SWE::n_variables);

            for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
                edge_internal.global_dof_indx[indx] = indx + global_dof_offset;
            }

            boundary.global_dof_indx = edge_internal.global_dof_indx;

            global_dof_offset += ndof_global * SWE::n_variables;
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

template <typename ProblemType>
void Problem::initialize_global_problem_parallel_finalize_pre_send(HDGDiscretization<ProblemType>& discretization,
                                                                   uint global_dof_offset) {
    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_offset](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint ndof_global = edge_int.edge_data.get_ndof();

        // Set indexes for global matrix construction
        for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
            edge_internal.global_dof_indx[indx] += global_dof_offset;
        }

        boundary_in.global_dof_indx = edge_internal.global_dof_indx;
        boundary_ex.global_dof_indx = edge_internal.global_dof_indx;
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        uint ndof_global = edge_bound.edge_data.get_ndof();

        // Set indexes for global matrix construction
        for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
            edge_internal.global_dof_indx[indx] += global_dof_offset;
        }

        boundary.global_dof_indx = edge_internal.global_dof_indx;
    });

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&global_dof_offset](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            // Set indexes for global matrix construction
            for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
                edge_internal.global_dof_indx[indx] += global_dof_offset;
            }

            boundary.global_dof_indx = edge_internal.global_dof_indx;

            std::vector<double> message(1);

            message[0] = (double)edge_internal.global_dof_indx[0];

            edge_dbound.boundary.boundary_condition.exchanger.SetToSendBuffer(CommTypes::global_dof_indx, message);
        }
    });
}

template <typename ProblemType, typename Communicator>
void Problem::initialize_global_problem_parallel_post_receive(HDGDiscretization<ProblemType>& discretization,
                                                              Communicator& communicator,
                                                              std::vector<uint>& global_dof_indx) {
    discretization.mesh_skeleton.CallForEachEdgeInterface([&global_dof_indx](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        edge_internal.delta_hat_global_con_flat.resize(edge_int.interface.data_in.get_nbound() +
                                                       edge_int.interface.data_ex.get_nbound() - 2);

        edge_internal.sol_offset = global_dof_indx.size();

        global_dof_indx.insert(
            global_dof_indx.end(), edge_internal.global_dof_indx.begin(), edge_internal.global_dof_indx.end());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&global_dof_indx](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        edge_internal.delta_hat_global_con_flat.resize(edge_bound.boundary.data.get_nbound() - 1);

        edge_internal.sol_offset = global_dof_indx.size();

        global_dof_indx.insert(
            global_dof_indx.end(), edge_internal.global_dof_indx.begin(), edge_internal.global_dof_indx.end());
    });

    discretization.mesh_skeleton.CallForEachEdgeDistributed([&global_dof_indx](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        edge_internal.delta_hat_global_con_flat.resize(edge_dbound.boundary.data.get_nbound() - 1);

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            std::vector<double> message(1);

            edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::global_dof_indx, message);

            uint global_dof_offset = (uint)message[0];

            edge_internal.global_dof_indx.resize(ndof_global * SWE::n_variables);

            // Set indexes for global matrix construction
            for (uint indx = 0; indx < ndof_global * SWE::n_variables; ++indx) {
                edge_internal.global_dof_indx[indx] = indx + global_dof_offset;
            }

            boundary.global_dof_indx = edge_internal.global_dof_indx;
        }

        edge_internal.sol_offset = global_dof_indx.size();

        global_dof_indx.insert(
            global_dof_indx.end(), edge_internal.global_dof_indx.begin(), edge_internal.global_dof_indx.end());

        edge_dbound.boundary.boundary_condition.ComputeInitTrace(
            edge_dbound, HybMatrix<double, SWE::n_variables>(), HybMatrix<double, SWE::n_variables>());
    });
}
}
}

#endif