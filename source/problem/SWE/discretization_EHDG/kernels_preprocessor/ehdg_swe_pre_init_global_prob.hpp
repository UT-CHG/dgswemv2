#ifndef EHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP
#define EHDG_SWE_PRE_INIT_GLOBAL_PROB_HPP

namespace SWE {
namespace EHDG {
template <typename ProblemType>
void Problem::initialize_global_problem_serial(HDGDiscretization<ProblemType>& discretization) {
    discretization.mesh_skeleton.CallForEachEdgeInterface([](auto& edge_int) {
        auto& edge_internal = edge_int.edge_data.edge_internal;

        auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
        auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

        uint gp_ex;
        for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
            gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

            edge_internal.aux_hat_at_gp(SWE::Auxiliaries::bath, gp) =
                (boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp) +
                 boundary_ex.aux_at_gp(SWE::Auxiliaries::bath, gp_ex)) /
                2.0;
        }

        set_constant(edge_int.edge_data.edge_state.q_hat, 0.0);

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof(),
                                              SWE::n_variables * edge_int.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_int.edge_data.get_ndof());
    });

    discretization.mesh_skeleton.CallForEachEdgeBoundary([](auto& edge_bound) {
        auto& edge_internal = edge_bound.edge_data.edge_internal;

        auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath) = row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        set_constant(edge_bound.edge_data.edge_state.q_hat, 0.0);

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_bound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_bound.edge_data.get_ndof());
    });
}

template <typename ProblemType>
void Problem::initialize_global_problem_parallel_pre_send(HDGDiscretization<ProblemType>& discretization) {
    Problem::initialize_global_problem_serial(discretization);

    discretization.mesh_skeleton.CallForEachEdgeDistributed([](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ngp = edge_dbound.boundary.data.get_ngp_boundary(edge_dbound.boundary.bound_id);

        // Construct message to exterior state
        std::vector<double> message;

        message.reserve(ngp);

        for (uint gp = 0; gp < ngp; ++gp) {
            message.push_back(boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));
        }

        edge_dbound.boundary.boundary_condition.exchanger.SetToSendBuffer(CommTypes::init_global_prob, message);

        set_constant(edge_dbound.edge_data.edge_state.q_hat, 0.0);

        // Initialize delta_hat_global and rhs_global containers
        edge_internal.delta_hat_global.resize(SWE::n_variables * edge_dbound.edge_data.get_ndof(),
                                              SWE::n_variables * edge_dbound.edge_data.get_ndof());
        edge_internal.rhs_global.resize(SWE::n_variables * edge_dbound.edge_data.get_ndof());
    });
}

template <typename ProblemType>
void Problem::initialize_global_problem_parallel_post_receive(HDGDiscretization<ProblemType>& discretization) {
    discretization.mesh_skeleton.CallForEachEdgeDistributed([](auto& edge_dbound) {
        auto& edge_internal = edge_dbound.edge_data.edge_internal;

        auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

        uint ngp = edge_dbound.boundary.data.get_ngp_boundary(edge_dbound.boundary.bound_id);

        std::vector<double> message(ngp);
        DynRowVector<double> bath_ex(ngp);

        edge_dbound.boundary.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::init_global_prob, message);

        for (uint gp = 0; gp < ngp; ++gp) {
            bath_ex[ngp - gp - 1] = message[gp];
        }

        row(edge_internal.aux_hat_at_gp, SWE::Auxiliaries::bath) =
            (row(boundary.aux_at_gp, SWE::Auxiliaries::bath) + bath_ex) / 2.0;
    });
}
}
}

#endif