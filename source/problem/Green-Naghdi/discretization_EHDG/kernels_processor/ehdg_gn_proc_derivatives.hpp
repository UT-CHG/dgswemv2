#ifndef EHDG_GN_PROC_DERIVATIVES_HPP
#define EHDG_GN_PROC_DERIVATIVES_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::serial_velocity_derivatives_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
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

        row(state.du, GN::FirstDerivatives::ux) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.u_at_gp, GlobalCoord::x));
        row(state.du, GN::FirstDerivatives::uy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.u_at_gp, GlobalCoord::x));
        row(state.du, GN::FirstDerivatives::vx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.u_at_gp, GlobalCoord::y));
        row(state.du, GN::FirstDerivatives::vy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.u_at_gp, GlobalCoord::y));
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

        /* IN State */
        row(boundary_in.u_hat_at_gp, GlobalCoord::x) =
            cwise_division(row(boundary_in.q_at_gp, GN::Variables::qx), row(boundary_in.aux_at_gp, GN::Auxiliaries::h));
        row(boundary_in.u_hat_at_gp, GlobalCoord::y) =
            cwise_division(row(boundary_in.q_at_gp, GN::Variables::qy), row(boundary_in.aux_at_gp, GN::Auxiliaries::h));

        /* EX State */
        row(boundary_ex.u_hat_at_gp, GlobalCoord::x) =
            cwise_division(row(boundary_ex.q_at_gp, GN::Variables::qx), row(boundary_ex.aux_at_gp, GN::Auxiliaries::h));
        row(boundary_ex.u_hat_at_gp, GlobalCoord::y) =
            cwise_division(row(boundary_ex.q_at_gp, GN::Variables::qy), row(boundary_ex.aux_at_gp, GN::Auxiliaries::h));

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            boundary_in.u_hat_at_gp(GlobalCoord::x, gp) =
                (boundary_in.u_hat_at_gp(GlobalCoord::x, gp) + boundary_ex.u_hat_at_gp(GlobalCoord::x, gp_ex)) / 2.0;

            boundary_in.u_hat_at_gp(GlobalCoord::y, gp) =
                (boundary_in.u_hat_at_gp(GlobalCoord::y, gp) + boundary_ex.u_hat_at_gp(GlobalCoord::y, gp_ex)) / 2.0;

            boundary_ex.u_hat_at_gp(GlobalCoord::x, gp_ex) = boundary_in.u_hat_at_gp(GlobalCoord::x, gp);

            boundary_ex.u_hat_at_gp(GlobalCoord::y, gp_ex) = boundary_in.u_hat_at_gp(GlobalCoord::y, gp);
        }

        /* IN State */
        row(state_in.du, GN::FirstDerivatives::ux) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.du, GN::FirstDerivatives::uy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.du, GN::FirstDerivatives::vx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.du, GN::FirstDerivatives::vy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_in, GlobalCoord::y)));

        /* EX State */
        row(state_ex.du, GN::FirstDerivatives::ux) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.du, GN::FirstDerivatives::uy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.du, GN::FirstDerivatives::vx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.du, GN::FirstDerivatives::vy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_ex, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.q_at_gp = bound.ComputeUgp(state.q);

        row(boundary.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary.q_at_gp, GN::Variables::ze) + row(boundary.aux_at_gp, GN::Auxiliaries::bath);

        row(boundary.u_hat_at_gp, GlobalCoord::x) =
            cwise_division(row(boundary.q_at_gp, GN::Variables::qx), row(boundary.aux_at_gp, GN::Auxiliaries::h));
        row(boundary.u_hat_at_gp, GlobalCoord::y) =
            cwise_division(row(boundary.q_at_gp, GN::Variables::qy), row(boundary.aux_at_gp, GN::Auxiliaries::h));

        row(state.du, GN::FirstDerivatives::ux) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::x), row(bound.surface_normal, GlobalCoord::x)));
        row(state.du, GN::FirstDerivatives::uy) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::x), row(bound.surface_normal, GlobalCoord::y)));
        row(state.du, GN::FirstDerivatives::vx) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::y), row(bound.surface_normal, GlobalCoord::x)));
        row(state.du, GN::FirstDerivatives::vy) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::y), row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        state.du = elt.ApplyMinv(state.du);
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state    = elt.data.state[stage];
        auto& internal = elt.data.internal;

        internal.du_at_gp = elt.ComputeUgp(state.du);

        row(state.ddu, GN::SecondDerivatives::uxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::FirstDerivatives::ux));
        row(state.ddu, GN::SecondDerivatives::uxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::FirstDerivatives::ux));
        row(state.ddu, GN::SecondDerivatives::uyx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::FirstDerivatives::uy));
        row(state.ddu, GN::SecondDerivatives::uyy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::FirstDerivatives::uy));

        row(state.ddu, GN::SecondDerivatives::vxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::FirstDerivatives::vx));
        row(state.ddu, GN::SecondDerivatives::vxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::FirstDerivatives::vx));
        row(state.ddu, GN::SecondDerivatives::vyx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::FirstDerivatives::vy));
        row(state.ddu, GN::SecondDerivatives::vyy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::FirstDerivatives::vy));
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

            boundary_in.du_hat_at_gp(GN::FirstDerivatives::ux, gp) =
                (boundary_in.du_hat_at_gp(GN::FirstDerivatives::ux, gp) +
                 boundary_ex.du_hat_at_gp(GN::FirstDerivatives::ux, gp_ex)) /
                2.0;

            boundary_in.du_hat_at_gp(GN::FirstDerivatives::uy, gp) =
                (boundary_in.du_hat_at_gp(GN::FirstDerivatives::uy, gp) +
                 boundary_ex.du_hat_at_gp(GN::FirstDerivatives::uy, gp_ex)) /
                2.0;

            boundary_in.du_hat_at_gp(GN::FirstDerivatives::vx, gp) =
                (boundary_in.du_hat_at_gp(GN::FirstDerivatives::vx, gp) +
                 boundary_ex.du_hat_at_gp(GN::FirstDerivatives::vx, gp_ex)) /
                2.0;

            boundary_in.du_hat_at_gp(GN::FirstDerivatives::vy, gp) =
                (boundary_in.du_hat_at_gp(GN::FirstDerivatives::vy, gp) +
                 boundary_ex.du_hat_at_gp(GN::FirstDerivatives::vy, gp_ex)) /
                2.0;

            boundary_ex.du_hat_at_gp(GN::FirstDerivatives::ux, gp_ex) =
                boundary_in.du_hat_at_gp(GN::FirstDerivatives::ux, gp);

            boundary_ex.du_hat_at_gp(GN::FirstDerivatives::uy, gp_ex) =
                boundary_in.du_hat_at_gp(GN::FirstDerivatives::uy, gp);

            boundary_ex.du_hat_at_gp(GN::FirstDerivatives::vx, gp_ex) =
                boundary_in.du_hat_at_gp(GN::FirstDerivatives::vx, gp);

            boundary_ex.du_hat_at_gp(GN::FirstDerivatives::vy, gp_ex) =
                boundary_in.du_hat_at_gp(GN::FirstDerivatives::vy, gp);
        }

        /* IN State */
        row(state_in.ddu, GN::SecondDerivatives::uxx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::ux), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::SecondDerivatives::uxy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::ux), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.ddu, GN::SecondDerivatives::uyx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::uy), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::SecondDerivatives::uyy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::uy), row(intface.surface_normal_in, GlobalCoord::y)));

        row(state_in.ddu, GN::SecondDerivatives::vxx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::vx), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::SecondDerivatives::vxy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::vx), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.ddu, GN::SecondDerivatives::vyx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::vy), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::SecondDerivatives::vyy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::FirstDerivatives::vy), row(intface.surface_normal_in, GlobalCoord::y)));
        /* EX State */
        row(state_ex.ddu, GN::SecondDerivatives::uxx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::ux), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::SecondDerivatives::uxy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::ux), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.ddu, GN::SecondDerivatives::uyx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::uy), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::SecondDerivatives::uyy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::uy), row(intface.surface_normal_ex, GlobalCoord::y)));

        row(state_ex.ddu, GN::SecondDerivatives::vxx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::vx), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::SecondDerivatives::vxy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::vx), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.ddu, GN::SecondDerivatives::vyx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::vy), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::SecondDerivatives::vyy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::FirstDerivatives::vy), row(intface.surface_normal_ex, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.du_hat_at_gp = bound.ComputeUgp(state.du);

        row(state.ddu, GN::SecondDerivatives::uxx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::ux), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::SecondDerivatives::uxy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::ux), row(bound.surface_normal, GlobalCoord::y)));
        row(state.ddu, GN::SecondDerivatives::uyx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::uy), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::SecondDerivatives::uyy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::uy), row(bound.surface_normal, GlobalCoord::y)));

        row(state.ddu, GN::SecondDerivatives::vxx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::vx), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::SecondDerivatives::vxy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::vx), row(bound.surface_normal, GlobalCoord::y)));
        row(state.ddu, GN::SecondDerivatives::vyx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::vy), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::SecondDerivatives::vyy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::FirstDerivatives::vy), row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        state.ddu = elt.ApplyMinv(state.ddu);
    });
}

void Problem::serial_bathymetry_derivatives_kernel(const RKStepper& stepper,
                                                   ProblemDiscretizationType& discretization) {}
}
}

#endif