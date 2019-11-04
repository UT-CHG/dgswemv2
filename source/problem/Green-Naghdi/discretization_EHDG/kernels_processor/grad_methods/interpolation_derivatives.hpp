#ifndef INTERPOLATION_DERIVATIVES_HPP
#define INTERPOLATION_DERIVATIVES_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void interpolate_dze(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
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
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
        std::vector<double> message(GN::n_dimensions + GN::n_du_terms + ngp);
        dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
        for (uint dim = 0; dim < GN::n_dimensions; ++dim) {
            derivative.dze_at_baryctr_neigh[dbound.bound_id][dim] = message[dim];
        }
        for (uint du = 0; du < GN::n_du_terms; ++du) {
            derivative.du_at_baryctr_neigh[dbound.bound_id][du] = message[GN::n_dimensions + du];
        }
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            boundary.aux_at_gp(SWE::Auxiliaries::h, gp) =
                (boundary.aux_at_gp(SWE::Auxiliaries::h, gp) + message[GN::n_dimensions + GN::n_du_terms + gp_ex]) /
                2.0;
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;

        for (uint bound = 0; bound < elt.data.get_nbound(); ++bound) {
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.dze_at_midpts, bound) =
                derivative.dze_at_baryctr +
                (derivative.dze_at_baryctr_neigh[element_1] - derivative.dze_at_baryctr) * derivative.a[bound][0] +
                (derivative.dze_at_baryctr_neigh[element_2] - derivative.dze_at_baryctr) * derivative.a[bound][1];
        }
        derivative.dze_lin = derivative.dze_at_midpts * derivative.T;
        state.dze          = elt.ProjectLinearToBasis(derivative.dze_lin);
    });
}

template <typename ProblemDiscretizationType>
void interpolate_du(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
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
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.du_at_midpts, bound) =
                derivative.du_at_baryctr +
                (derivative.du_at_baryctr_neigh[element_1] - derivative.du_at_baryctr) * derivative.a[bound][0] +
                (derivative.du_at_baryctr_neigh[element_2] - derivative.du_at_baryctr) * derivative.a[bound][1];
        }
        derivative.du_lin = derivative.du_at_midpts * derivative.T;
        state.du          = elt.ProjectLinearToBasis(derivative.du_lin);
        internal.du_at_gp = elt.ComputeUgp(state.du);
    });
}

template <typename ProblemDiscretizationType>
void interpolate_ddu(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
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
            const uint element_1 = derivative.a_elem[2 * bound];
            const uint element_2 = derivative.a_elem[2 * bound + 1];
            column(derivative.ddu_at_midpts, bound) =
                derivative.ddu_at_baryctr +
                (derivative.ddu_at_baryctr_neigh[element_1] - derivative.ddu_at_baryctr) * derivative.a[bound][0] +
                (derivative.ddu_at_baryctr_neigh[element_2] - derivative.ddu_at_baryctr) * derivative.a[bound][1];
        }
        derivative.ddu_lin = derivative.ddu_at_midpts * derivative.T;
        state.ddu          = elt.ProjectLinearToBasis(derivative.ddu_lin);
        internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
    });
}
}
}

#endif
