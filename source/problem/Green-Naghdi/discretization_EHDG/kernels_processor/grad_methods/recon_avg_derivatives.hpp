#ifndef RECON_AVG_DERIVATIVES_HPP
#define RECON_AVG_DERIVATIVES_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_dze(ProblemDiscretizationType& discretization,
                     ProblemGlobalDataType& global_data,
                     const ESSPRKStepper& stepper) {
    set_constant(global_data.derivatives_at_node, 0.0);

    discretization.mesh.CallForEachElement([&global_data](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            subvector(
                global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_dimensions, GN::n_dimensions) +=
                derivative.dze_at_baryctr;
        }
    });

    discretization.mesh.CallForEachElement([&global_data, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            column(derivative.dze_lin, node) = subvector(global_data.derivatives_at_node,
                                                         derivative.local_nodeID[node] * GN::n_dimensions,
                                                         GN::n_dimensions) /
                                               derivative.node_mult[node];
        }
        state.dze = elt.ProjectLinearToBasis(derivative.dze_lin);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_du(ProblemDiscretizationType& discretization,
                    ProblemGlobalDataType& global_data,
                    const ESSPRKStepper& stepper) {
    set_constant(global_data.derivatives_at_node, 0.0);

    discretization.mesh.CallForEachElement([&global_data](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            subvector(
                global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_du_terms, GN::n_du_terms) +=
                derivative.du_at_baryctr;
        }
    });

    discretization.mesh.CallForEachElement([&global_data, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            column(derivative.du_lin, node) =
                subvector(
                    global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_du_terms, GN::n_du_terms) /
                derivative.node_mult[node];
        }
        state.du          = elt.ProjectLinearToBasis(derivative.du_lin);
        internal.du_at_gp = elt.ComputeUgp(state.du);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_ddu(ProblemDiscretizationType& discretization,
                     ProblemGlobalDataType& global_data,
                     const ESSPRKStepper& stepper) {
    set_constant(global_data.derivatives_at_node, 0.0);

    discretization.mesh.CallForEachElement([&global_data](auto& elt) {
        auto& derivative = elt.data.derivative;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            subvector(
                global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_ddu_terms, GN::n_ddu_terms) +=
                derivative.ddu_at_baryctr;
        }
    });

    discretization.mesh.CallForEachElement([&global_data, &stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;
        for (uint node = 0; node < elt.GetNodeID().size(); ++node) {
            column(derivative.ddu_lin, node) =
                subvector(
                    global_data.derivatives_at_node, derivative.local_nodeID[node] * GN::n_ddu_terms, GN::n_ddu_terms) /
                derivative.node_mult[node];
        }
        state.ddu          = elt.ProjectLinearToBasis(derivative.ddu_lin);
        internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
    });
}
}
}

#endif
