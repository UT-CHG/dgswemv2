#ifndef EHDG_GN_PRE_DBATH_SERIAL_HPP
#define EHDG_GN_PRE_DBATH_SERIAL_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
void Problem::serial_bathymetry_derivatives_kernel(ProblemDiscretizationType& discretization) {
    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        row(state.dbath, GlobalCoord::x) = -elt.IntegrationDPhi(GlobalCoord::x, row(internal.aux_at_gp, GN::Auxiliaries::bath));
        row(state.dbath, GlobalCoord::y) = -elt.IntegrationDPhi(GlobalCoord::y, row(internal.aux_at_gp, GN::Auxiliaries::bath));
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& state_in    = intface.data_in.state[0];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[0];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            boundary_in.bath_hat_at_gp[gp] =
                (boundary_in.aux_at_gp(GN::Auxiliaries::bath, gp) + boundary_ex.aux_at_gp(GN::Auxiliaries::bath, gp_ex)) / 2.0;

            boundary_ex.bath_hat_at_gp[gp_ex] = boundary_in.bath_hat_at_gp[gp];
        }

        /* IN State */
        row(state_in.dbath, GlobalCoord::x) += intface.IntegrationPhiIN(cwise_multiplication(
            boundary_in.bath_hat_at_gp, row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.dbath, GlobalCoord::y) += intface.IntegrationPhiIN(cwise_multiplication(
            boundary_in.bath_hat_at_gp, row(intface.surface_normal_in, GlobalCoord::y)));

        /* EX State */
        row(state_ex.dbath, GlobalCoord::x) += intface.IntegrationPhiEX(cwise_multiplication(
            boundary_ex.bath_hat_at_gp, row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.dbath, GlobalCoord::y) += intface.IntegrationPhiEX(cwise_multiplication(
            boundary_ex.bath_hat_at_gp, row(intface.surface_normal_ex, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state    = bound.data.state[0];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.bath_hat_at_gp = row(boundary.aux_at_gp, GN::Auxiliaries::bath);

        row(state.dbath, GlobalCoord::x) += bound.IntegrationPhi(
            cwise_multiplication(boundary.bath_hat_at_gp, row(bound.surface_normal, GlobalCoord::x)));
        row(state.dbath, GlobalCoord::y) += bound.IntegrationPhi(
            cwise_multiplication(boundary.bath_hat_at_gp, row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];

        state.dbath = elt.ApplyMinv(state.dbath);
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state    = elt.data.state[0];
        auto& internal = elt.data.internal;

        internal.dbath_at_gp = elt.ComputeUgp(state.dbath);

        row(state.ddbath, GN::DDBath::bxx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.dbath_at_gp, GlobalCoord::x));
        row(state.ddbath, GN::DDBath::bxy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.dbath_at_gp, GlobalCoord::x));
        row(state.ddbath, GN::DDBath::byx) =
            -elt.IntegrationDPhi(GlobalCoord::x, row(internal.dbath_at_gp, GlobalCoord::y));
        row(state.ddbath, GN::DDBath::byy) =
            -elt.IntegrationDPhi(GlobalCoord::y, row(internal.dbath_at_gp, GlobalCoord::y));
    });

    discretization.mesh.CallForEachInterface([](auto& intface) {
        auto& state_in    = intface.data_in.state[0];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[0];
        auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

        boundary_in.dbath_at_gp = intface.ComputeUgpIN(state_in.dbath);
        boundary_ex.dbath_at_gp = intface.ComputeUgpEX(state_ex.dbath);

        uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
        uint gp_ex;
        for (uint gp = 0; gp < ngp; ++gp) {
            gp_ex = ngp - gp - 1;

            boundary_in.dbath_hat_at_gp(GlobalCoord::x, gp) =
                (boundary_in.dbath_at_gp(GlobalCoord::x, gp) + boundary_ex.dbath_at_gp(GlobalCoord::x, gp_ex)) / 2.0;
            boundary_in.dbath_hat_at_gp(GlobalCoord::y, gp) =
                (boundary_in.dbath_at_gp(GlobalCoord::y, gp) + boundary_ex.dbath_at_gp(GlobalCoord::y, gp_ex)) / 2.0;

            boundary_ex.dbath_hat_at_gp(GlobalCoord::x, gp_ex) = boundary_in.dbath_hat_at_gp(GlobalCoord::x, gp);
            boundary_ex.dbath_hat_at_gp(GlobalCoord::y, gp_ex) = boundary_in.dbath_hat_at_gp(GlobalCoord::y, gp);
        }

        /* IN State */
        row(state_in.ddbath, GN::DDBath::bxx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.dbath_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddbath, GN::DDBath::bxy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.dbath_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_in, GlobalCoord::y)));
        row(state_in.ddbath, GN::DDBath::byx) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.dbath_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_in, GlobalCoord::x)));
        row(state_in.ddbath, GN::DDBath::byy) += intface.IntegrationPhiIN(cwise_multiplication(
            row(boundary_in.dbath_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_in, GlobalCoord::y)));

        /* EX State */
        row(state_ex.ddbath, GN::DDBath::bxx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.dbath_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddbath, GN::DDBath::bxy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.dbath_hat_at_gp, GlobalCoord::x), row(intface.surface_normal_ex, GlobalCoord::y)));
        row(state_ex.ddbath, GN::DDBath::byx) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.dbath_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_ex, GlobalCoord::x)));
        row(state_ex.ddbath, GN::DDBath::byy) += intface.IntegrationPhiEX(cwise_multiplication(
            row(boundary_ex.dbath_hat_at_gp, GlobalCoord::y), row(intface.surface_normal_ex, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachBoundary([](auto& bound) {
        auto& state    = bound.data.state[0];
        auto& boundary = bound.data.boundary[bound.bound_id];

        boundary.dbath_at_gp = bound.ComputeUgp(state.dbath);

        boundary.dbath_hat_at_gp = boundary.dbath_at_gp;

        row(state.ddbath, GN::DDBath::bxx) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.dbath_at_gp, GlobalCoord::x), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddbath, GN::DDBath::bxy) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.dbath_at_gp, GlobalCoord::x), row(bound.surface_normal, GlobalCoord::y)));
        row(state.ddbath, GN::DDBath::byx) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.dbath_at_gp, GlobalCoord::y), row(bound.surface_normal, GlobalCoord::x)));
        row(state.ddbath, GN::DDBath::byy) += bound.IntegrationPhi(
            cwise_multiplication(row(boundary.dbath_at_gp, GlobalCoord::y), row(bound.surface_normal, GlobalCoord::y)));
    });

    discretization.mesh.CallForEachElement([](auto& elt) {
        auto& state = elt.data.state[0];

        state.ddbath = elt.ApplyMinv(state.ddbath);
    });

    /*discretization.mesh.CallForEachElement([](auto& elt) {


        auto& state    = elt.data.state[0];
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

    discretization.mesh.CallForEachInterface([](auto& intface) {


        auto& state_in    = intface.data_in.state[0];
        auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

        auto& state_ex    = intface.data_ex.state[0];
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
        }*/

        /* IN State */
        /*row(state_in.ddu, GN::DDU::uxx) += intface.IntegrationPhiIN(cwise_multiplication(
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
            row(boundary_in.du_hat_at_gp, GN::DU::vy), row(intface.surface_normal_in, GlobalCoord::y)));*/
        /* EX State */
        /*row(state_ex.ddu, GN::DDU::uxx) += intface.IntegrationPhiEX(cwise_multiplication(
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

    discretization.mesh.CallForEachBoundary([](auto& bound) {


        auto& state    = bound.data.state[0];
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

    discretization.mesh.CallForEachElement([](auto& elt) {


        auto& state = elt.data.state[0];

        state.ddu = elt.ApplyMinv(state.ddu);
    });*/
}
}
}

#endif