#ifndef IHDG_GN_PROC_DERIVATIVES_HPP
#define IHDG_GN_PROC_DERIVATIVES_HPP

#include "general_definitions.hpp"

namespace GN {
namespace IHDG {
void Problem::serial_derivatives_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.q_at_gp = elt.ComputeUgp(state.q);

        row(internal.aux_at_gp, GN::Auxiliaries::h) =
            row(internal.q_at_gp, GN::Variables::ze) + row(internal.aux_at_gp, GN::Auxiliaries::bath);

        row(internal.u_at_gp, GlobalCoord::x) =
            cwise_division(row(internal.q_at_gp, GN::Variables::qx), row(internal.aux_at_gp, GN::Auxiliaries::h));
        row(internal.u_at_gp, GlobalCoord::y) =
            cwise_division(row(internal.q_at_gp, GN::Variables::qy), row(internal.aux_at_gp, GN::Auxiliaries::h));

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state.dze, dir) = -elt.IntegrationDPhi(dir, row(internal.q_at_gp, GN::Variables::ze));
        }

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.du, GN::n_dimensions * u + dir) = -elt.IntegrationDPhi(dir, row(internal.u_at_gp, u));
            }
        }
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage = stepper.GetStage();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);

        row(boundary_in.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary_in.q_at_gp, GN::Variables::ze) + row(boundary_in.aux_at_gp, GN::Auxiliaries::bath);

        row(boundary_ex.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary_ex.q_at_gp, GN::Variables::ze) + row(boundary_ex.aux_at_gp, GN::Auxiliaries::bath);

        row(boundary_in.u_hat_at_gp, GlobalCoord::x) =
            cwise_division(row(boundary_in.q_at_gp, GN::Variables::qx), row(boundary_in.aux_at_gp, GN::Auxiliaries::h));
        row(boundary_in.u_hat_at_gp, GlobalCoord::y) =
            cwise_division(row(boundary_in.q_at_gp, GN::Variables::qy), row(boundary_in.aux_at_gp, GN::Auxiliaries::h));

        row(boundary_ex.u_hat_at_gp, GlobalCoord::x) =
            cwise_division(row(boundary_ex.q_at_gp, GN::Variables::qx), row(boundary_ex.aux_at_gp, GN::Auxiliaries::h));
        row(boundary_ex.u_hat_at_gp, GlobalCoord::y) =
            cwise_division(row(boundary_ex.q_at_gp, GN::Variables::qy), row(boundary_ex.aux_at_gp, GN::Auxiliaries::h));

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            boundary_in.ze_hat_at_gp[gp] =
                (boundary_in.q_at_gp(GN::Variables::ze, gp) + boundary_ex.q_at_gp(GN::Variables::ze, gp_ex)) / 2.0;

            boundary_ex.ze_hat_at_gp[gp_ex] = boundary_in.ze_hat_at_gp[gp];

            for (uint u = 0; u < GN::n_dimensions; ++u) {
                boundary_in.u_hat_at_gp(u, gp) =
                    (boundary_in.u_hat_at_gp(u, gp) + boundary_ex.u_hat_at_gp(u, gp_ex)) / 2.0;

                boundary_ex.u_hat_at_gp(u, gp_ex) = boundary_in.u_hat_at_gp(u, gp);
            }
        }

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state_in.dze, dir) += intface.IntegrationPhiIN(
                cwise_multiplication(boundary_in.ze_hat_at_gp, row(intface.surface_normal_in, dir)));

            row(state_ex.dze, dir) += intface.IntegrationPhiEX(
                cwise_multiplication(boundary_ex.ze_hat_at_gp, row(intface.surface_normal_ex, dir)));
        }

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.du, GN::n_dimensions * u + dir) += intface.IntegrationPhiIN(
                    cwise_multiplication(row(boundary_in.u_hat_at_gp, u), row(intface.surface_normal_in, dir)));

                row(state_ex.du, GN::n_dimensions * u + dir) += intface.IntegrationPhiEX(
                    cwise_multiplication(row(boundary_ex.u_hat_at_gp, u), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.q_at_gp = bound.ComputeUgp(state.q);

        row(boundary.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary.q_at_gp, GN::Variables::ze) + row(boundary.aux_at_gp, GN::Auxiliaries::bath);

        boundary.ze_hat_at_gp = row(boundary.q_at_gp, GN::Variables::ze);

        row(boundary.u_hat_at_gp, GlobalCoord::x) =
            cwise_division(row(boundary.q_at_gp, GN::Variables::qx), row(boundary.aux_at_gp, GN::Auxiliaries::h));
        row(boundary.u_hat_at_gp, GlobalCoord::y) =
            cwise_division(row(boundary.q_at_gp, GN::Variables::qy), row(boundary.aux_at_gp, GN::Auxiliaries::h));

        for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
            row(state.dze, dir) +=
                bound.IntegrationPhi(cwise_multiplication(boundary.ze_hat_at_gp, row(bound.surface_normal, dir)));
        }

        for (uint u = 0; u < GN::n_dimensions; ++u) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.du, GN::n_dimensions * u + dir) += bound.IntegrationPhi(
                    cwise_multiplication(row(boundary.u_hat_at_gp, u), row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        state.dze = elt.ApplyMinv(state.dze);
        state.du  = elt.ApplyMinv(state.du);
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.du_at_gp = elt.ComputeUgp(state.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.ddu, GN::n_dimensions * du + dir) = -elt.IntegrationDPhi(dir, row(internal.du_at_gp, du));
            }
        }
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage = stepper.GetStage();

        auto& state_in    = intface.data_in.state[stage];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[stage];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.du_hat_at_gp = intface.ComputeUgpIN(state_in.du);
        boundary_ex.du_hat_at_gp = intface.ComputeUgpEX(state_ex.du);

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            for (uint du = 0; du < GN::n_du_terms; ++du) {
                boundary_in.du_hat_at_gp(du, gp) =
                    (boundary_in.du_hat_at_gp(du, gp) + boundary_ex.du_hat_at_gp(du, gp_ex)) / 2.0;

                boundary_ex.du_hat_at_gp(du, gp_ex) = boundary_in.du_hat_at_gp(du, gp);
            }
        }

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.ddu, GN::n_dimensions * du + dir) += intface.IntegrationPhiIN(
                    cwise_multiplication(row(boundary_in.du_hat_at_gp, du), row(intface.surface_normal_in, dir)));

                row(state_ex.ddu, GN::n_dimensions * du + dir) += intface.IntegrationPhiEX(
                    cwise_multiplication(row(boundary_ex.du_hat_at_gp, du), row(intface.surface_normal_ex, dir)));
            }
        }
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.du_hat_at_gp = bound.ComputeUgp(state.du);

        for (uint du = 0; du < GN::n_du_terms; ++du) {
            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.ddu, GN::n_dimensions * du + dir) += bound.IntegrationPhi(
                    cwise_multiplication(row(boundary.du_hat_at_gp, du), row(bound.surface_normal, dir)));
            }
        }
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        state.ddu = elt.ApplyMinv(state.ddu);

        internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
    });
}
}
}

#endif