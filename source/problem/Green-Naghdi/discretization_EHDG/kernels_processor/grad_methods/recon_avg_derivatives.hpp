#ifndef RECON_AVG_DERIVATIVES_HPP
#define RECON_AVG_DERIVATIVES_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_dze(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data, const ESSPRKStepper& stepper) {
    set_constant(global_data.derivatives_at_node, 0.0);

    discretization.mesh.CallForEachElement([&global_data](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            subvector(global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_dimensions, GN::n_dimensions) += derivative.dze_at_baryctr;
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&global_data](auto& dbound) {
        auto& derivative = dbound.data.derivative;

        std::vector<double> message(GN::n_dimensions + GN::n_du_terms);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        for (uint dim = 0; dim < GN::n_dimensions; ++dim) {
            derivative.dze_at_baryctr_neigh[dbound.bound_id][dim] = message[dim];
        }
        for (uint du = 0; du < GN::n_du_terms; ++du) {
            derivative.du_at_baryctr_neigh[dbound.bound_id][du] = message[GN::n_dimensions + du];
        }

        for (uint i = 0; i < 2; ++i) {
            subvector(global_data.derivatives_at_node,
                      derivative.local_nodeID[(dbound.bound_id + i + 1) % 3] *  GN::n_dimensions, GN::n_dimensions) +=
                        derivative.dze_at_baryctr_neigh[dbound.bound_id];
        }
    });

    discretization.mesh.CallForEachElement([&global_data, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;

        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            column(derivative.dze_lin, node) =
                    subvector(global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_dimensions, GN::n_dimensions) /
                            derivative.node_mult[node];
        }
        state.dze = elt.ProjectLinearToBasis(derivative.dze_lin);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_du(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data, const ESSPRKStepper& stepper) {
    set_constant(global_data.derivatives_at_node, 0.0);

    discretization.mesh.CallForEachElement([&global_data](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            subvector(global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_du_terms, GN::n_du_terms) += derivative.du_at_baryctr;
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&global_data](auto& dbound) {
        auto& derivative = dbound.data.derivative;

        for (uint i = 0; i < 2; ++i) {
            subvector(global_data.derivatives_at_node,
                      derivative.local_nodeID[(dbound.bound_id + i + 1) % 3] *  GN::n_du_terms, GN::n_du_terms) +=
                    derivative.du_at_baryctr_neigh[dbound.bound_id];
        }
    });

    discretization.mesh.CallForEachElement([&global_data, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            column(derivative.du_lin, node) =
                    subvector(global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_du_terms, GN::n_du_terms) /
                            derivative.node_mult[node];
        }
        state.du          = elt.ProjectLinearToBasis(derivative.du_lin);
        internal.du_at_gp = elt.ComputeUgp(state.du);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_ddu(ProblemDiscretizationType& discretization, ProblemGlobalDataType& global_data, const ESSPRKStepper& stepper) {
    set_constant(global_data.derivatives_at_node, 0.0);

    discretization.mesh.CallForEachElement([&global_data](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            subvector(global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_ddu_terms, GN::n_ddu_terms) += derivative.ddu_at_baryctr;
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&global_data](auto& dbound) {
        auto& derivative = dbound.data.derivative;

        std::vector<double> message(n_ddu_terms);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        for (uint ddu = 0; ddu < GN::n_ddu_terms; ++ddu) {
            derivative.ddu_at_baryctr_neigh[dbound.bound_id][ddu] = message[ddu];
        }

        for (uint i = 0; i < 2; ++i) {
            subvector(global_data.derivatives_at_node,
                      derivative.local_nodeID[(dbound.bound_id + i + 1) % 3] *  GN::n_ddu_terms, GN::n_ddu_terms) +=
                    derivative.ddu_at_baryctr_neigh[dbound.bound_id];
        }
    });

    discretization.mesh.CallForEachElement([&global_data, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            column(derivative.ddu_lin, node) =
                    subvector(global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_ddu_terms, GN::n_ddu_terms) /
                            derivative.node_mult[node];
        }
        state.ddu          = elt.ProjectLinearToBasis(derivative.ddu_lin);
        internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
    });
}
}
}

#endif
