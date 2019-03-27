#ifndef IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define IHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace IHDG {
template <typename ProblemType>
void Problem::initialize_global_problem_serial(HDGDiscretization<ProblemType>& discretization,
                                               const typename ProblemType::ProblemStepperType& stepper,
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

        uint gp_ex;
        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

            edge_internal.aux_hat_at_gp(SWE::Auxiliaries::bath, gp) =
                (boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp) +
                 boundary_ex.aux_at_gp(SWE::Auxiliaries::bath, gp_ex)) /
                2.0;
        }

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

    discretization.mesh_skeleton.CallForEachEdgeBoundary([&stepper, &global_dof_offset](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath) = row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

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

        edge_bound.boundary.boundary_condition.ComputeInitTrace(stepper, edge_bound);
    });
}

template <typename ProblemType>
void Problem::initialize_global_problem_parallel_pre_send(HDGDiscretization<ProblemType>& discretization,
                                                          const typename ProblemType::ProblemStepperType& stepper,
                                                          uint& global_dof_offset) {
    Problem::initialize_global_problem_serial(discretization, stepper, global_dof_offset);

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
        auto& dbound = edge_dbound.boundary;

        auto& state    = dbound.data.state[0];
        auto& boundary = dbound.data.boundary[dbound.bound_id];

        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        boundary.q_at_gp = dbound.ComputeUgp(state.q);

        // Construct message to exterior state
        std::vector<double> message;

        message.reserve(1 + (SWE::n_variables + 1) * dbound.data.get_ngp_boundary(dbound.bound_id));

        for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
            for (uint var = 0; var < SWE::n_variables; ++var) {
                message.push_back(boundary.q_at_gp(var, gp));
            }
        }

        for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
            message.push_back(boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));
        }

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

            message.push_back((double)edge_internal.global_dof_indx[0]);
        }

        edge_dbound.boundary.boundary_condition.exchanger.SetToSendBuffer(CommTypes::init_global_prob, message);
    });
}

template <typename ProblemType>
void Problem::initialize_global_problem_parallel_post_receive(HDGDiscretization<ProblemType>& discretization,
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
        auto& dbound = edge_dbound.boundary;

        auto& boundary = dbound.data.boundary[dbound.bound_id];

        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        uint ngp = edge_dbound.boundary.data.get_ngp_boundary(edge_dbound.boundary.bound_id);

        std::vector<double> message(1 + (SWE::n_variables + 1) * ngp);
        HybMatrix<double, SWE::n_variables> q_ex(SWE::n_variables, ngp);
        DynRowVector<double> bath_ex(ngp);

        edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::init_global_prob, message);

        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            for (uint var = 0; var < SWE::n_variables; ++var) {
                q_ex(var, gp_ex) = message[SWE::n_variables * gp + var];
            }
        }

        for (uint gp = 0; gp < ngp; ++gp) {
            bath_ex(ngp - gp - 1) = message[SWE::n_variables * ngp + gp];
        }

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath) =
            (row(boundary.aux_at_gp, SWE::Auxiliaries::bath) + bath_ex) / 2.0;

        edge_internal.delta_hat_global_con_flat.resize(edge_dbound.boundary.data.get_nbound() - 1);

        uint ndof_global = edge_dbound.edge_data.get_ndof();

        uint locality_in = edge_dbound.boundary.boundary_condition.exchanger.locality_in;
        uint submesh_in  = edge_dbound.boundary.boundary_condition.exchanger.submesh_in;
        uint locality_ex = edge_dbound.boundary.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = edge_dbound.boundary.boundary_condition.exchanger.submesh_ex;

        if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            uint global_dof_offset = (uint)message.back();

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

        edge_dbound.boundary.boundary_condition.ComputeInitTrace(edge_dbound, q_ex);
    });
}
}
}

#endif