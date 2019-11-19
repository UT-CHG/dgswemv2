#ifndef RECON_LS_DERIVATIVES_HPP
#define RECON_LS_DERIVATIVES_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_dze(ProblemDiscretizationType& discretization,
                     ProblemGlobalDataType& global_data,
                     const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& derivative_in                                     = intface.data_in.derivative;
        auto& derivative_ex                                     = intface.data_ex.derivative;
        derivative_in.dze_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.dze_at_baryctr;
        derivative_ex.dze_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.dze_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& derivative                                = bound.data.derivative;
        derivative.dze_at_baryctr_neigh[bound.bound_id] = derivative.dze_at_baryctr;
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        std::vector<double> message(GN::n_dimensions + GN::n_du_terms);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        for (uint dim = 0; dim < GN::n_dimensions; ++dim) {
            derivative.dze_at_baryctr_neigh[dbound.bound_id][dim] = message[dim];
        }
        for (uint du = 0; du < GN::n_du_terms; ++du) {
            derivative.du_at_baryctr_neigh[dbound.bound_id][du] = message[GN::n_dimensions + du];
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.dze_at_midpts, bound) =
                derivative.dze_at_baryctr_neigh[bound] - derivative.dze_at_baryctr;
        }
        derivative.dze_at_midpts *= transpose(derivative.P) * derivative.dX_transpose;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.dze_at_midpts, bound) += derivative.dze_at_baryctr;
        }
        derivative.dze_lin = derivative.dze_at_midpts * derivative.T;
        state.dze          = elt.ProjectLinearToBasis(derivative.dze_lin);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_du(ProblemDiscretizationType& discretization,
                    ProblemGlobalDataType& global_data,
                    const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& derivative_in                                    = intface.data_in.derivative;
        auto& derivative_ex                                    = intface.data_ex.derivative;
        derivative_in.du_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.du_at_baryctr;
        derivative_ex.du_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.du_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& derivative                               = bound.data.derivative;
        derivative.du_at_baryctr_neigh[bound.bound_id] = derivative.du_at_baryctr;
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.du_at_midpts, bound) = derivative.du_at_baryctr_neigh[bound] - derivative.du_at_baryctr;
        }
        derivative.du_at_midpts *= transpose(derivative.P) * derivative.dX_transpose;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.du_at_midpts, bound) += derivative.du_at_baryctr;
        }
        derivative.du_lin = derivative.du_at_midpts * derivative.T;
        state.du          = elt.ProjectLinearToBasis(derivative.du_lin);
        internal.du_at_gp = elt.ComputeUgp(state.du);
    });
}

template <typename ProblemDiscretizationType, typename ProblemGlobalDataType>
void reconstruct_ddu(ProblemDiscretizationType& discretization,
                     ProblemGlobalDataType& global_data,
                     const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& derivative_in                                     = intface.data_in.derivative;
        auto& derivative_ex                                     = intface.data_ex.derivative;
        derivative_in.ddu_at_baryctr_neigh[intface.bound_id_in] = derivative_ex.ddu_at_baryctr;
        derivative_ex.ddu_at_baryctr_neigh[intface.bound_id_ex] = derivative_in.ddu_at_baryctr;
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& derivative                                = bound.data.derivative;
        derivative.ddu_at_baryctr_neigh[bound.bound_id] = derivative.ddu_at_baryctr;
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        std::vector<double> message(n_ddu_terms);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        for (uint ddu = 0; ddu < GN::n_ddu_terms; ++ddu) {
            derivative.ddu_at_baryctr_neigh[dbound.bound_id][ddu] = message[ddu];
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.ddu_at_midpts, bound) =
                derivative.ddu_at_baryctr_neigh[bound] - derivative.ddu_at_baryctr;
        }
        derivative.ddu_at_midpts *= transpose(derivative.P) * derivative.dX_transpose;
        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            column(derivative.ddu_at_midpts, bound) += derivative.ddu_at_baryctr;
        }
        derivative.ddu_lin = derivative.ddu_at_midpts * derivative.T;
        state.ddu          = elt.ProjectLinearToBasis(derivative.ddu_lin);
        internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
    });
}
}
}

#endif
