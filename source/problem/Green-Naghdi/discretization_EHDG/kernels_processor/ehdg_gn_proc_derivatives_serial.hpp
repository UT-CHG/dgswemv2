#ifndef EHDG_GN_PROC_DERIVATIVES_SERIAL_HPP
#define EHDG_GN_PROC_DERIVATIVES_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dze_avg(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_du_avg(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_ddu_avg(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);

template <typename ProblemDiscretizationType>
void compute_dze(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_du(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_ddu(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);

void Problem::compute_derivatives_serial(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    compute_dze_avg(discretization, stepper);
    compute_dze(discretization, stepper);

    compute_du_avg(discretization, stepper);
    compute_du(discretization, stepper);

    compute_ddu_avg(discretization, stepper);
    compute_ddu(discretization, stepper);
}

template <typename ProblemDiscretizationType>
void compute_dze_avg(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);
        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        set_constant(derivative.dze_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage    = stepper.GetStage();
        auto& state_in      = intface.data_in.state[stage];
        auto& state_ex      = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);
        row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative_in.dze_at_baryctr, dir) +=
                1.0 / derivative_in.area *
                intface.IntegrationIN(
                    vec_cw_mult(row(boundary_in.q_at_gp, SWE::Variables::ze), row(intface.surface_normal_in, dir)));
            row(derivative_ex.dze_at_baryctr, dir) +=
                1.0 / derivative_ex.area *
                intface.IntegrationEX(
                    vec_cw_mult(row(boundary_ex.q_at_gp, SWE::Variables::ze), row(intface.surface_normal_ex, dir)));
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        boundary.q_at_gp = bound.ComputeUgp(state.q);
        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative.dze_at_baryctr, dir) +=
                1.0 / derivative.area *
                bound.Integration(
                    vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(bound.surface_normal, dir)));
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        const uint stage = stepper.GetStage();
        auto& state      = dbound.data.state[stage];
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        boundary.q_at_gp = dbound.ComputeUgp(state.q);
        row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(derivative.dze_at_baryctr, dir) +=
                1.0 / derivative.area *
                dbound.Integration(
                    vec_cw_mult(row(boundary.q_at_gp, SWE::Variables::ze), row(dbound.surface_normal, dir)));
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_dze(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
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
void compute_du_avg(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& derivative = elt.data.derivative;
        auto& internal   = elt.data.internal;

        row(internal.u_at_gp, GlobalCoord::x) =
            vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
        row(internal.u_at_gp, GlobalCoord::y) =
            vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

        set_constant(derivative.du_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in   = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex   = intface.data_ex.boundary[intface.bound_id_ex];

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative_in.du_at_baryctr, GN::n_dimensions * u + dir) +=
                    1.0 / derivative_in.area *
                    intface.IntegrationIN(vec_cw_mult(vec_cw_div(row(boundary_in.q_at_gp, SWE::Variables::qx + u),
                                                                 row(boundary_in.aux_at_gp, SWE::Auxiliaries::h)),
                                                      row(intface.surface_normal_in, dir)));
                row(derivative_ex.du_at_baryctr, GN::n_dimensions * u + dir) +=
                    1.0 / derivative_ex.area *
                    intface.IntegrationEX(vec_cw_mult(vec_cw_div(row(boundary_ex.q_at_gp, SWE::Variables::qx + u),
                                                                 row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h)),
                                                      row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.du_at_baryctr, GN::n_dimensions * u + dir) +=
                    1.0 / derivative.area *
                    bound.Integration(vec_cw_mult(vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx + u),
                                                             row(boundary.aux_at_gp, SWE::Auxiliaries::h)),
                                                  row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        auto& derivative = dbound.data.derivative;
        auto& boundary   = dbound.data.boundary[dbound.bound_id];

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.du_at_baryctr, GN::n_dimensions * u + dir) +=
                    1.0 / derivative.area *
                    dbound.Integration(vec_cw_mult(vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx + u),
                                                              row(boundary.aux_at_gp, SWE::Auxiliaries::h)),
                                                   row(dbound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_du(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
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
void compute_ddu_avg(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& derivative = elt.data.derivative;
        set_constant(derivative.ddu_at_baryctr, 0);
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage    = stepper.GetStage();
        auto& state_in      = intface.data_in.state[stage];
        auto& state_ex      = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;

        const auto du_at_gp_in = intface.ComputeUgpIN(state_in.du);
        const auto du_at_gp_ex = intface.ComputeUgpEX(state_ex.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative_in.ddu_at_baryctr, GN::n_dimensions * du + dir) +=
                    1.0 / derivative_in.area *
                    intface.IntegrationIN(vec_cw_mult(row(du_at_gp_in, du), row(intface.surface_normal_in, dir)));
                row(derivative_ex.ddu_at_baryctr, GN::n_dimensions * du + dir) +=
                    1.0 / derivative_ex.area *
                    intface.IntegrationEX(vec_cw_mult(row(du_at_gp_ex, du), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;

        const auto du_at_gp = bound.ComputeUgp(state.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.ddu_at_baryctr, GN::n_dimensions * du + dir) +=
                    1.0 / derivative.area *
                    bound.Integration(vec_cw_mult(row(du_at_gp, du), row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
        const uint stage = stepper.GetStage();
        auto& state      = dbound.data.state[stage];
        auto& derivative = dbound.data.derivative;

        const auto du_at_gp = dbound.ComputeUgp(state.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(derivative.ddu_at_baryctr, GN::n_dimensions * du + dir) +=
                    1.0 / derivative.area *
                    dbound.Integration(vec_cw_mult(row(du_at_gp, du), row(dbound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_ddu(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
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
