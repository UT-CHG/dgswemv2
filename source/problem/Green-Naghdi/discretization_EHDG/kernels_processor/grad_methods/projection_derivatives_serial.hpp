#ifndef PROJECTION_DERIVATIVES_SERIAL_HPP
#define PROJECTION_DERIVATIVES_SERIAL_HPP

namespace GN {
namespace EHDG {
template <typename ProblemDiscretizationType>
void compute_dze_rhs(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_du_rhs(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);
template <typename ProblemDiscretizationType>
void compute_ddu_rhs(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper);

void Problem::compute_derivatives_serial(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    compute_dze_rhs(discretization, stepper);
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        state.dze        = elt.ApplyMinv(state.dze);
    });

    compute_du_rhs(discretization, stepper);
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        state.du         = elt.ApplyMinv(state.du);
    });

    compute_ddu_rhs(discretization, stepper);
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage   = stepper.GetStage();
        auto& state        = elt.data.state[stage];
        auto& internal     = elt.data.internal;
        state.ddu          = elt.ApplyMinv(state.ddu);
        internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
    });
}

template <typename ProblemDiscretizationType>
void compute_dze_rhs(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& internal   = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);
        row(internal.aux_at_gp, SWE::Auxiliaries::h) =
            row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state.dze, dir) = -elt.IntegrationDPhi(dir, row(internal.q_at_gp, SWE::Variables::ze));
        }
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage  = stepper.GetStage();
        auto& state_in    = intface.data_in.state[stage];
        auto& state_ex    = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);
        row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);
        row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
            row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);
        derivative_in.ze_hat_at_gp[intface.bound_id_in] = row(boundary_in.q_at_gp, SWE::Variables::ze);
        derivative_ex.ze_hat_at_gp[intface.bound_id_ex] = row(boundary_ex.q_at_gp, SWE::Variables::ze);

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            derivative_in.ze_hat_at_gp[intface.bound_id_in][gp] =
                (derivative_in.ze_hat_at_gp[intface.bound_id_in][gp] + derivative_ex.ze_hat_at_gp[intface.bound_id_ex][gp_ex]) / 2.0;
            derivative_ex.ze_hat_at_gp[intface.bound_id_ex][gp_ex] = derivative_in.ze_hat_at_gp[intface.bound_id_in][gp];
        }

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state_in.dze, dir) +=
                intface.IntegrationPhiIN(vec_cw_mult(derivative_in.ze_hat_at_gp[intface.bound_id_in], row(intface.surface_normal_in, dir)));
            row(state_ex.dze, dir) +=
                intface.IntegrationPhiEX(vec_cw_mult(derivative_ex.ze_hat_at_gp[intface.bound_id_ex], row(intface.surface_normal_ex, dir)));
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
        derivative.ze_hat_at_gp[bound.bound_id] = row(boundary.q_at_gp, SWE::Variables::ze);

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state.dze, dir) +=
                bound.IntegrationPhi(vec_cw_mult(derivative.ze_hat_at_gp[bound.bound_id], row(bound.surface_normal, dir)));
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_du_rhs(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& internal   = elt.data.internal;

        row(internal.u_at_gp, GlobalCoord::x) =
            vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
        row(internal.u_at_gp, GlobalCoord::y) =
            vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.du, GN::n_dimensions * u + dir) = -elt.IntegrationDPhi(dir, row(internal.u_at_gp, u));
            }
        }
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage  = stepper.GetStage();
        auto& state_in    = intface.data_in.state[stage];
        auto& state_ex    = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        row(derivative_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::x) =
                vec_cw_div(row(boundary_in.q_at_gp, SWE::Variables::qx), row(boundary_in.aux_at_gp, SWE::Auxiliaries::h));
        row(derivative_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::y) =
                vec_cw_div(row(boundary_in.q_at_gp, SWE::Variables::qy), row(boundary_in.aux_at_gp, SWE::Auxiliaries::h));

        row(derivative_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::x) =
                vec_cw_div(row(boundary_ex.q_at_gp, SWE::Variables::qx), row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h));
        row(derivative_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::y) =
                vec_cw_div(row(boundary_ex.q_at_gp, SWE::Variables::qy), row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h));

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint u = 0; u < GN::n_dimensions; ++u) {
                derivative_in.u_hat_at_gp[intface.bound_id_in](u, gp) =
                    (derivative_in.u_hat_at_gp[intface.bound_id_in](u, gp) +
                     derivative_ex.u_hat_at_gp[intface.bound_id_ex](u, gp_ex)) / 2.0;
                derivative_ex.u_hat_at_gp[intface.bound_id_ex](u, gp_ex) = derivative_in.u_hat_at_gp[intface.bound_id_in](u, gp);
            }
        }

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.du, GN::n_dimensions * u + dir) += intface.IntegrationPhiIN(
                    vec_cw_mult(row(derivative_in.u_hat_at_gp[intface.bound_id_in], u), row(intface.surface_normal_in, dir)));
                row(state_ex.du, GN::n_dimensions * u + dir) += intface.IntegrationPhiEX(
                    vec_cw_mult(row(derivative_ex.u_hat_at_gp[intface.bound_id_ex], u), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;
        auto& boundary   = bound.data.boundary[bound.bound_id];

        row(derivative.u_hat_at_gp[bound.bound_id], GlobalCoord::x) =
            vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
        row(derivative.u_hat_at_gp[bound.bound_id], GlobalCoord::y) =
            vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.du, GN::n_dimensions * u + dir) +=
                    bound.IntegrationPhi(vec_cw_mult(row(derivative.u_hat_at_gp[bound.bound_id], u), row(bound.surface_normal, dir)));
            }
        }
    });
}

template <typename ProblemDiscretizationType>
void compute_ddu_rhs(ProblemDiscretizationType& discretization, const ESSPRKStepper& stepper) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();
        auto& state      = elt.data.state[stage];
        auto& internal   = elt.data.internal;

        internal.du_at_gp = elt.ComputeUgp(state.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.ddu, GN::n_dimensions * du + dir) = -elt.IntegrationDPhi(dir, row(internal.du_at_gp, du));
            }
        }
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage  = stepper.GetStage();
        auto& state_in    = intface.data_in.state[stage];
        auto& state_ex    = intface.data_ex.state[stage];
        auto& derivative_in = intface.data_in.derivative;
        auto& derivative_ex = intface.data_ex.derivative;

        derivative_in.du_hat_at_gp[intface.bound_id_in] = intface.ComputeUgpIN(state_in.du);
        derivative_ex.du_hat_at_gp[intface.bound_id_ex] = intface.ComputeUgpEX(state_ex.du);

        const uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        for (uint gp = 0; gp < ngp; ++gp) {
            const uint gp_ex = ngp - gp - 1;
            for (uint du = 0; du < GN::n_du_terms; ++du) {
                derivative_in.du_hat_at_gp[intface.bound_id_in](du, gp) =
                    (derivative_in.du_hat_at_gp[intface.bound_id_in](du, gp) +
                     derivative_ex.du_hat_at_gp[intface.bound_id_ex](du, gp_ex)) / 2.0;
                derivative_ex.du_hat_at_gp[intface.bound_id_ex](du, gp_ex) = derivative_in.du_hat_at_gp[intface.bound_id_in](du, gp);
            }
        }

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.ddu, GN::n_dimensions * du + dir) += intface.IntegrationPhiIN(
                    vec_cw_mult(row(derivative_in.du_hat_at_gp[intface.bound_id_in], du), row(intface.surface_normal_in, dir)));
                row(state_ex.ddu, GN::n_dimensions * du + dir) += intface.IntegrationPhiEX(
                    vec_cw_mult(row(derivative_ex.du_hat_at_gp[intface.bound_id_ex], du), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();
        auto& state      = bound.data.state[stage];
        auto& derivative = bound.data.derivative;

        derivative.du_hat_at_gp[bound.bound_id] = bound.ComputeUgp(state.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.ddu, GN::n_dimensions * du + dir) +=
                    bound.IntegrationPhi(vec_cw_mult(row(derivative.du_hat_at_gp[bound.bound_id], du), row(bound.surface_normal, dir)));
            }
        }
    });
}
}
}

#endif