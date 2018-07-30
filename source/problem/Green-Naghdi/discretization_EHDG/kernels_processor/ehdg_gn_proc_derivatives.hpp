#ifndef EHDG_GN_PROC_DERIVATIVES_HPP
#define EHDG_GN_PROC_DERIVATIVES_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::serial_derivatives_kernel(const RKStepper& stepper, ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        const uint stage = stepper.GetStage();

        auto& state     = elt.data.state[stage];
        auto& internal  = elt.data.internal;
        auto& disp_corr = elt.data.disp_corr;

        internal.q_at_gp = elt.ComputeUgp(state.q);

        row(internal.aux_at_gp, GN::Auxiliaries::h) =
            row(internal.q_at_gp, GN::Variables::ze) + row(internal.aux_at_gp, GN::Auxiliaries::bath);

        row(disp_corr.u_at_gp, GlobalCoord::x) =
            cwise_division(row(internal.q_at_gp, GN::Variables::qx), row(internal.aux_at_gp, GN::Auxiliaries::h));
        row(disp_corr.u_at_gp, GlobalCoord::y) =
            cwise_division(row(internal.q_at_gp, GN::Variables::qy), row(internal.aux_at_gp, GN::Auxiliaries::h));

        row(disp_corr.du, GN::FirstDerivatives::ux) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(disp_corr.u_at_gp, GlobalCoord::x));
        row(disp_corr.du, GN::FirstDerivatives::uy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(disp_corr.u_at_gp, GlobalCoord::x));
        row(disp_corr.du, GN::FirstDerivatives::vx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(disp_corr.u_at_gp, GlobalCoord::y));
        row(disp_corr.du, GN::FirstDerivatives::vy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(disp_corr.u_at_gp, GlobalCoord::y));
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        const uint stage = stepper.GetStage();

        auto& state_in     = intface.data_in.state[stage];
        auto& boundary_in  = intface.data_in.boundary[intface.bound_id_in];
        auto& disp_corr_in = intface.data_in.disp_corr;

        auto& state_ex     = intface.data_ex.state[stage];
        auto& boundary_ex  = intface.data_ex.boundary[intface.bound_id_ex];
        auto& disp_corr_ex = intface.data_ex.disp_corr;

        boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
        boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);

        row(boundary_in.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary_in.q_at_gp, GN::Variables::ze) + row(boundary_in.aux_at_gp, GN::Auxiliaries::bath);

        row(boundary_ex.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary_ex.q_at_gp, GN::Variables::ze) + row(boundary_ex.aux_at_gp, GN::Auxiliaries::bath);

        /* IN State */
        row(disp_corr_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::x) =
            cwise_division(row(boundary_in.q_at_gp, GN::Variables::qx), row(boundary_in.aux_at_gp, GN::Auxiliaries::h));
        row(disp_corr_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::y) =
            cwise_division(row(boundary_in.q_at_gp, GN::Variables::qy), row(boundary_in.aux_at_gp, GN::Auxiliaries::h));

        /* EX State */
        row(disp_corr_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::x) =
            cwise_division(row(boundary_ex.q_at_gp, GN::Variables::qx), row(boundary_ex.aux_at_gp, GN::Auxiliaries::h));
        row(disp_corr_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::y) =
            cwise_division(row(boundary_ex.q_at_gp, GN::Variables::qy), row(boundary_ex.aux_at_gp, GN::Auxiliaries::h));

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            disp_corr_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::x, gp) =
                (disp_corr_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::x, gp) +
                 disp_corr_ex.u_hat_at_gp[intface.bound_id_ex](GlobalCoord::x, gp_ex)) /
                2.0;

            disp_corr_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::y, gp) =
                (disp_corr_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::y, gp) +
                 disp_corr_ex.u_hat_at_gp[intface.bound_id_ex](GlobalCoord::y, gp_ex)) /
                2.0;

            disp_corr_ex.u_hat_at_gp[intface.bound_id_ex](GlobalCoord::x, gp_ex) =
                disp_corr_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::x, gp);

            disp_corr_ex.u_hat_at_gp[intface.bound_id_ex](GlobalCoord::y, gp_ex) =
                disp_corr_in.u_hat_at_gp[intface.bound_id_in](GlobalCoord::y, gp);
        }

        /* IN State */
        row(disp_corr_in.du, GN::FirstDerivatives::ux) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::x),
                                 row(intface.surface_normal_in, GlobalCoord::x)));
        row(disp_corr_in.du, GN::FirstDerivatives::uy) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::x),
                                 row(intface.surface_normal_in, GlobalCoord::y)));
        row(disp_corr_in.du, GN::FirstDerivatives::vx) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::y),
                                 row(intface.surface_normal_in, GlobalCoord::x)));
        row(disp_corr_in.du, GN::FirstDerivatives::vy) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.u_hat_at_gp[intface.bound_id_in], GlobalCoord::y),
                                 row(intface.surface_normal_in, GlobalCoord::y)));

        /* EX State */
        row(disp_corr_ex.du, GN::FirstDerivatives::ux) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::x),
                                 row(intface.surface_normal_ex, GlobalCoord::x)));
        row(disp_corr_ex.du, GN::FirstDerivatives::uy) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::x),
                                 row(intface.surface_normal_ex, GlobalCoord::y)));
        row(disp_corr_ex.du, GN::FirstDerivatives::vx) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::y),
                                 row(intface.surface_normal_ex, GlobalCoord::x)));
        row(disp_corr_ex.du, GN::FirstDerivatives::vy) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.u_hat_at_gp[intface.bound_id_ex], GlobalCoord::y),
                                 row(intface.surface_normal_ex, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        const uint stage = stepper.GetStage();

        auto& state     = bound.data.state[stage];
        auto& boundary  = bound.data.boundary[bound.bound_id];
        auto& disp_corr = bound.data.disp_corr;

        boundary.q_at_gp = bound.ComputeUgp(state.q);

        row(boundary.aux_at_gp, GN::Auxiliaries::h) =
            row(boundary.q_at_gp, GN::Variables::ze) + row(boundary.aux_at_gp, GN::Auxiliaries::bath);

        row(disp_corr.u_hat_at_gp[bound.bound_id], GlobalCoord::x) =
            cwise_division(row(boundary.q_at_gp, GN::Variables::qx), row(boundary.aux_at_gp, GN::Auxiliaries::h));
        row(disp_corr.u_hat_at_gp[bound.bound_id], GlobalCoord::y) =
            cwise_division(row(boundary.q_at_gp, GN::Variables::qy), row(boundary.aux_at_gp, GN::Auxiliaries::h));

        row(disp_corr.du, GN::FirstDerivatives::ux) += bound.IntegrationPhi(cwise_multiplication(
            row(disp_corr.u_hat_at_gp[bound.bound_id], GlobalCoord::x), row(bound.surface_normal, GlobalCoord::x)));
        row(disp_corr.du, GN::FirstDerivatives::uy) += bound.IntegrationPhi(cwise_multiplication(
            row(disp_corr.u_hat_at_gp[bound.bound_id], GlobalCoord::x), row(bound.surface_normal, GlobalCoord::y)));
        row(disp_corr.du, GN::FirstDerivatives::vx) += bound.IntegrationPhi(cwise_multiplication(
            row(disp_corr.u_hat_at_gp[bound.bound_id], GlobalCoord::y), row(bound.surface_normal, GlobalCoord::x)));
        row(disp_corr.du, GN::FirstDerivatives::vy) += bound.IntegrationPhi(cwise_multiplication(
            row(disp_corr.u_hat_at_gp[bound.bound_id], GlobalCoord::y), row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& disp_corr = elt.data.disp_corr;

        disp_corr.du = elt.ApplyMinv(disp_corr.du);
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& disp_corr = elt.data.disp_corr;

        disp_corr.du_at_gp = elt.ComputeUgp(disp_corr.du);

        row(disp_corr.ddu, GN::SecondDerivatives::uxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(disp_corr.du_at_gp, GN::FirstDerivatives::ux));
        row(disp_corr.ddu, GN::SecondDerivatives::uxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(disp_corr.du_at_gp, GN::FirstDerivatives::ux));
        row(disp_corr.ddu, GN::SecondDerivatives::uyx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(disp_corr.du_at_gp, GN::FirstDerivatives::uy));
        row(disp_corr.ddu, GN::SecondDerivatives::uyy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(disp_corr.du_at_gp, GN::FirstDerivatives::uy));

        row(disp_corr.ddu, GN::SecondDerivatives::vxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(disp_corr.du_at_gp, GN::FirstDerivatives::vx));
        row(disp_corr.ddu, GN::SecondDerivatives::vxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(disp_corr.du_at_gp, GN::FirstDerivatives::vx));
        row(disp_corr.ddu, GN::SecondDerivatives::vyx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(disp_corr.du_at_gp, GN::FirstDerivatives::vy));
        row(disp_corr.ddu, GN::SecondDerivatives::vyy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(disp_corr.du_at_gp, GN::FirstDerivatives::vy));
    });

    discretization.mesh.CallForEachInterface([&stepper](auto& intface) {
        auto& disp_corr_in = intface.data_in.disp_corr;
        auto& disp_corr_ex = intface.data_ex.disp_corr;

        disp_corr_in.du_hat_at_gp[intface.bound_id_in] = intface.ComputeUgpIN(disp_corr_in.du);
        disp_corr_ex.du_hat_at_gp[intface.bound_id_ex] = intface.ComputeUgpEX(disp_corr_ex.du);

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::ux, gp) =
                (disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::ux, gp) +
                 disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::ux, gp_ex)) /
                2.0;

            disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::uy, gp) =
                (disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::uy, gp) +
                 disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::uy, gp_ex)) /
                2.0;

            disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::vx, gp) =
                (disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::vx, gp) +
                 disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::vx, gp_ex)) /
                2.0;

            disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::vy, gp) =
                (disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::vy, gp) +
                 disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::vy, gp_ex)) /
                2.0;

            disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::ux, gp_ex) =
                disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::ux, gp);

            disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::uy, gp_ex) =
                disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::uy, gp);

            disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::vx, gp_ex) =
                disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::vx, gp);

            disp_corr_ex.du_hat_at_gp[intface.bound_id_ex](GN::FirstDerivatives::vy, gp_ex) =
                disp_corr_in.du_hat_at_gp[intface.bound_id_in](GN::FirstDerivatives::vy, gp);
        }

        /* IN State */
        row(disp_corr_in.ddu, GN::SecondDerivatives::uxx) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::ux),
                                 row(intface.surface_normal_in, GlobalCoord::x)));
        row(disp_corr_in.ddu, GN::SecondDerivatives::uxy) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::ux),
                                 row(intface.surface_normal_in, GlobalCoord::y)));
        row(disp_corr_in.ddu, GN::SecondDerivatives::uyx) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::uy),
                                 row(intface.surface_normal_in, GlobalCoord::x)));
        row(disp_corr_in.ddu, GN::SecondDerivatives::uyy) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::uy),
                                 row(intface.surface_normal_in, GlobalCoord::y)));

        row(disp_corr_in.ddu, GN::SecondDerivatives::vxx) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::vx),
                                 row(intface.surface_normal_in, GlobalCoord::x)));
        row(disp_corr_in.ddu, GN::SecondDerivatives::vxy) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::vx),
                                 row(intface.surface_normal_in, GlobalCoord::y)));
        row(disp_corr_in.ddu, GN::SecondDerivatives::vyx) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::vy),
                                 row(intface.surface_normal_in, GlobalCoord::x)));
        row(disp_corr_in.ddu, GN::SecondDerivatives::vyy) += intface.IntegrationPhiIN(
            cwise_multiplication(row(disp_corr_in.du_hat_at_gp[intface.bound_id_in], GN::FirstDerivatives::vy),
                                 row(intface.surface_normal_in, GlobalCoord::y)));
        /* EX State */
        row(disp_corr_ex.ddu, GN::SecondDerivatives::uxx) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::ux),
                                 row(intface.surface_normal_ex, GlobalCoord::x)));
        row(disp_corr_ex.ddu, GN::SecondDerivatives::uxy) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::ux),
                                 row(intface.surface_normal_ex, GlobalCoord::y)));
        row(disp_corr_ex.ddu, GN::SecondDerivatives::uyx) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::uy),
                                 row(intface.surface_normal_ex, GlobalCoord::x)));
        row(disp_corr_ex.ddu, GN::SecondDerivatives::uyy) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::uy),
                                 row(intface.surface_normal_ex, GlobalCoord::y)));

        row(disp_corr_ex.ddu, GN::SecondDerivatives::vxx) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::vx),
                                 row(intface.surface_normal_ex, GlobalCoord::x)));
        row(disp_corr_ex.ddu, GN::SecondDerivatives::vxy) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::vx),
                                 row(intface.surface_normal_ex, GlobalCoord::y)));
        row(disp_corr_ex.ddu, GN::SecondDerivatives::vyx) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::vy),
                                 row(intface.surface_normal_ex, GlobalCoord::x)));
        row(disp_corr_ex.ddu, GN::SecondDerivatives::vyy) += intface.IntegrationPhiEX(
            cwise_multiplication(row(disp_corr_ex.du_hat_at_gp[intface.bound_id_ex], GN::FirstDerivatives::vy),
                                 row(intface.surface_normal_ex, GlobalCoord::y)));                                 
    });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) {
        auto& disp_corr = bound.data.disp_corr;

        disp_corr.du_hat_at_gp[bound.bound_id] = bound.ComputeUgp(disp_corr.du);

        row(disp_corr.ddu, GN::SecondDerivatives::uxx) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::ux),
                                 row(bound.surface_normal, GlobalCoord::x)));
        row(disp_corr.ddu, GN::SecondDerivatives::uxy) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::ux),
                                 row(bound.surface_normal, GlobalCoord::y)));
        row(disp_corr.ddu, GN::SecondDerivatives::uyx) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::uy),
                                 row(bound.surface_normal, GlobalCoord::x)));
        row(disp_corr.ddu, GN::SecondDerivatives::uyy) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::uy),
                                 row(bound.surface_normal, GlobalCoord::y)));

        row(disp_corr.ddu, GN::SecondDerivatives::vxx) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::vx),
                                 row(bound.surface_normal, GlobalCoord::x)));
        row(disp_corr.ddu, GN::SecondDerivatives::vxy) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::vx),
                                 row(bound.surface_normal, GlobalCoord::y)));
        row(disp_corr.ddu, GN::SecondDerivatives::vyx) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::vy),
                                 row(bound.surface_normal, GlobalCoord::x)));
        row(disp_corr.ddu, GN::SecondDerivatives::vyy) += bound.IntegrationPhi(
            cwise_multiplication(row(disp_corr.du_hat_at_gp[bound.bound_id], GN::FirstDerivatives::vy),
                                 row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) {
        auto& disp_corr = elt.data.disp_corr;

        disp_corr.ddu = elt.ApplyMinv(disp_corr.ddu);
    });
}
}
}

#endif