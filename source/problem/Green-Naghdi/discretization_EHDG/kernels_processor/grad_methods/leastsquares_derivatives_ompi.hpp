#ifndef LEASTSQUARES_DERIVATIVES_OMPI_HPP
#define LEASTSQUARES_DERIVATIVES_OMPI_HPP

#include "leastsquares_derivatives_serial.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                       const ESSPRKStepper& stepper,
                                       const uint begin_sim_id,
                                       const uint end_sim_id) {
    /* First derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            const uint stage = stepper.GetStage();
            auto& state      = elt.data.state[stage];
            auto& derivative = elt.data.derivative;
            auto& internal   = elt.data.internal;

            internal.q_at_gp = elt.ComputeUgp(state.q);
            row(internal.aux_at_gp, SWE::Auxiliaries::h) =
                    row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

            derivative.ze_lin        = elt.ProjectBasisToLinear(row(state.q, SWE::Variables::ze));
            derivative.ze_at_baryctr = elt.ComputeLinearUbaryctr(derivative.ze_lin);
            derivative.ze_at_midpts  = elt.ComputeLinearUmidpts(derivative.ze_lin);

            row(internal.u_at_gp, GlobalCoord::x) =
                    vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
            row(internal.u_at_gp, GlobalCoord::y) =
                    vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

            derivative.u_lin        = elt.ProjectBasisToLinear(elt.L2Projection(internal.u_at_gp));
            derivative.u_at_baryctr = elt.ComputeLinearUbaryctr(derivative.u_lin);
            derivative.u_at_midpts  = elt.ComputeLinearUmidpts(derivative.u_lin);
        });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(1 + GN::n_dimensions);
            message[0] = derivative.ze_at_baryctr;
            message[1] = derivative.u_at_baryctr[GlobalCoord::x];
            message[2] = derivative.u_at_baryctr[GlobalCoord::y];
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives, stepper.GetTimestamp());

        compute_dze_ls(sim_units[su_id]->discretization, stepper);
        compute_du_ls(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            auto& boundary   = dbound.data.boundary[dbound.bound_id];

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_dimensions + GN::n_du_terms + ngp);
            for (uint dim = 0; dim < GN::n_dimensions; ++dim) {
                message[dim] = derivative.dze_at_baryctr[dim];
            }
            for (uint du = 0; du < GN::n_du_terms; ++du) {
                message[GN::n_dimensions + du] = derivative.du_at_baryctr[du];
            }
            for (uint gp = 0; gp < ngp; ++gp) {
                message[GN::n_dimensions + GN::n_du_terms + gp] = boundary.aux_at_gp(SWE::Auxiliaries::h, gp);
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives, stepper.GetTimestamp());

        interpolate_dze(sim_units[su_id]->discretization, stepper);
        interpolate_du(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }
    /* First derivatives end */

    /* Second derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            auto& derivative = elt.data.derivative;
            derivative.du_at_baryctr = elt.ComputeLinearUbaryctr(derivative.du_lin);
            derivative.du_at_midpts  = elt.ComputeLinearUmidpts(derivative.du_lin);
        });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(GN::n_du_terms);
            for (uint du = 0; du < GN::n_du_terms; ++du) {
                message[du] = derivative.du_at_baryctr[du];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives, stepper.GetTimestamp());

        compute_ddu_ls(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(GN::n_ddu_terms);
            for (uint ddu = 0; ddu < GN::n_ddu_terms; ++ddu) {
                message[ddu] = derivative.ddu_at_baryctr[ddu];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives, stepper.GetTimestamp());

        interpolate_ddu(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }
    /* Second derivatives end */
}
}
}

#endif