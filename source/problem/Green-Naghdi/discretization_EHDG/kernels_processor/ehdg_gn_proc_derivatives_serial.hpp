#ifndef EHDG_GN_PROC_DERIVATIVES_HPP
#define EHDG_GN_PROC_DERIVATIVES_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
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

        row(state.dze, GlobalCoord::x) = -elt.IntegrationDPhi(GlobalCoord::x, row(internal.q_at_gp, GN::Variables::ze));
        row(state.dze, GlobalCoord::y) = -elt.IntegrationDPhi(GlobalCoord::y, row(internal.q_at_gp, GN::Variables::ze));

        row(state.du, GN::DU::ux) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.u_at_gp, GlobalCoord::x));
        row(state.du, GN::DU::uy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.u_at_gp, GlobalCoord::x));
        row(state.du, GN::DU::vx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.u_at_gp, GlobalCoord::y));
        row(state.du, GN::DU::vy) =
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

            boundary_in.ze_hat_at_gp[gp] =
                (boundary_in.q_at_gp(GN::Variables::ze, gp) + boundary_ex.q_at_gp(GN::Variables::ze, gp_ex)) / 2.0;

            boundary_in.u_hat_at_gp(GlobalCoord::x, gp) =
                (boundary_in.u_hat_at_gp(GlobalCoord::x, gp) + boundary_ex.u_hat_at_gp(GlobalCoord::x, gp_ex)) / 2.0;
            boundary_in.u_hat_at_gp(GlobalCoord::y, gp) =
                (boundary_in.u_hat_at_gp(GlobalCoord::y, gp) + boundary_ex.u_hat_at_gp(GlobalCoord::y, gp_ex)) / 2.0;

            boundary_ex.ze_hat_at_gp[gp_ex] = boundary_in.ze_hat_at_gp[gp];

            boundary_ex.u_hat_at_gp(GlobalCoord::x, gp_ex) = boundary_in.u_hat_at_gp(GlobalCoord::x, gp);
            boundary_ex.u_hat_at_gp(GlobalCoord::y, gp_ex) = boundary_in.u_hat_at_gp(GlobalCoord::y, gp);
        }

        /* IN State */
        row(state_in.dze, GlobalCoord::x) += intface.IntegrationPhiIN(cwise_multiplication(
            boundary_in.ze_hat_at_gp, row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.dze, GlobalCoord::y) += intface.IntegrationPhiIN(cwise_multiplication(
            boundary_in.ze_hat_at_gp, row(intface.surface_normal_in, GlobalCoord::y)));

        row(state_in.du, GN::DU::ux) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.du, GN::DU::uy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.du, GN::DU::vx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.du, GN::DU::vy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_in, GlobalCoord::y)));

        /* EX State */
        row(state_ex.dze, GlobalCoord::x) += intface.IntegrationPhiEX(cwise_multiplication(
            boundary_ex.ze_hat_at_gp, row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.dze, GlobalCoord::y) += intface.IntegrationPhiEX(cwise_multiplication(
            boundary_ex.ze_hat_at_gp, row(intface.surface_normal_ex, GlobalCoord::y)));

        row(state_ex.du, GN::DU::ux) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.du, GN::DU::uy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.du, GN::DU::vx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.du, GN::DU::vy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.u_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_ex, GlobalCoord::y)));
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

        row(state.dze, GlobalCoord::x) += bound.IntegrationPhi(
            cwise_multiplication(boundary.ze_hat_at_gp, row(bound.surface_normal, GlobalCoord::x)));
        row(state.dze, GlobalCoord::y) += bound.IntegrationPhi(
            cwise_multiplication(boundary.ze_hat_at_gp, row(bound.surface_normal, GlobalCoord::y)));

        row(state.du, GN::DU::ux) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::x), row(bound.surface_normal, GlobalCoord::x)));
        row(state.du, GN::DU::uy) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::x), row(bound.surface_normal, GlobalCoord::y)));
        row(state.du, GN::DU::vx) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::y), row(bound.surface_normal, GlobalCoord::x)));
        row(state.du, GN::DU::vy) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.u_hat_at_gp, GlobalCoord::y), row(bound.surface_normal, GlobalCoord::y)));
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

        row(state.ddu, GN::DDU::uxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::DU::ux));
        row(state.ddu, GN::DDU::uxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::DU::ux));
        row(state.ddu, GN::DDU::uyx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::DU::uy));
        row(state.ddu, GN::DDU::uyy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::DU::uy));

        row(state.ddu, GN::DDU::vxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::DU::vx));
        row(state.ddu, GN::DDU::vxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::DU::vx));
        row(state.ddu, GN::DDU::vyx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.du_at_gp, GN::DU::vy));
        row(state.ddu, GN::DDU::vyy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.du_at_gp, GN::DU::vy));
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

            boundary_in.du_hat_at_gp(GN::DU::ux, gp) =
                (boundary_in.du_hat_at_gp(GN::DU::ux, gp) +
                 boundary_ex.du_hat_at_gp(GN::DU::ux, gp_ex)) /
                2.0;

            boundary_in.du_hat_at_gp(GN::DU::uy, gp) =
                (boundary_in.du_hat_at_gp(GN::DU::uy, gp) +
                 boundary_ex.du_hat_at_gp(GN::DU::uy, gp_ex)) /
                2.0;

            boundary_in.du_hat_at_gp(GN::DU::vx, gp) =
                (boundary_in.du_hat_at_gp(GN::DU::vx, gp) +
                 boundary_ex.du_hat_at_gp(GN::DU::vx, gp_ex)) /
                2.0;

            boundary_in.du_hat_at_gp(GN::DU::vy, gp) =
                (boundary_in.du_hat_at_gp(GN::DU::vy, gp) +
                 boundary_ex.du_hat_at_gp(GN::DU::vy, gp_ex)) /
                2.0;

            boundary_ex.du_hat_at_gp(GN::DU::ux, gp_ex) =
                boundary_in.du_hat_at_gp(GN::DU::ux, gp);

            boundary_ex.du_hat_at_gp(GN::DU::uy, gp_ex) =
                boundary_in.du_hat_at_gp(GN::DU::uy, gp);

            boundary_ex.du_hat_at_gp(GN::DU::vx, gp_ex) =
                boundary_in.du_hat_at_gp(GN::DU::vx, gp);

            boundary_ex.du_hat_at_gp(GN::DU::vy, gp_ex) =
                boundary_in.du_hat_at_gp(GN::DU::vy, gp);
        }

        /* IN State */
        row(state_in.ddu, GN::DDU::uxx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::ux), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::DDU::uxy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::ux), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.ddu, GN::DDU::uyx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::uy), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::DDU::uyy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::uy), row(intface.surface_normal_in, GlobalCoord::y)));

        row(state_in.ddu, GN::DDU::vxx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::vx), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::DDU::vxy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::vx), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.ddu, GN::DDU::vyx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::vy), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddu, GN::DDU::vyy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.du_hat_at_gp, GN::DU::vy), row(intface.surface_normal_in, GlobalCoord::y)));
        /* EX State */
        row(state_ex.ddu, GN::DDU::uxx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::ux), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::DDU::uxy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::ux), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.ddu, GN::DDU::uyx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::uy), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::DDU::uyy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::uy), row(intface.surface_normal_ex, GlobalCoord::y)));

        row(state_ex.ddu, GN::DDU::vxx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::vx), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::DDU::vxy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::vx), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.ddu, GN::DDU::vyx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::vy), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddu, GN::DDU::vyy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.du_hat_at_gp, GN::DU::vy), row(intface.surface_normal_ex, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();

        auto& state    = bound.data.state[stage];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.du_hat_at_gp = bound.ComputeUgp(state.du);

        row(state.ddu, GN::DDU::uxx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::ux), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::DDU::uxy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::ux), row(bound.surface_normal, GlobalCoord::y)));
        row(state.ddu, GN::DDU::uyx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::uy), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::DDU::uyy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::uy), row(bound.surface_normal, GlobalCoord::y)));

        row(state.ddu, GN::DDU::vxx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::vx), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::DDU::vxy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::vx), row(bound.surface_normal, GlobalCoord::y)));
        row(state.ddu, GN::DDU::vyx) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::vy), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddu, GN::DDU::vyy) += bound.IntegrationPhi(cwise_multiplication(
            row(boundary.du_hat_at_gp, GN::DU::vy), row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state = elt.data.state[stage];

        state.ddu = elt.ApplyMinv(state.ddu);
    });
}
}
}

#endif