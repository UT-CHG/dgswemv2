#ifndef EHDG_GN_PROC_DERIVATIVES_OMPI_HPP
#define EHDG_GN_PROC_DERIVATIVES_OMPI_HPP

#include "ehdg_gn_proc_derivatives_serial.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                       const ESSPRKStepper& stepper,
                                       const uint begin_sim_id,
                                       const uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        compute_dze_avg(sim_units[su_id]->discretization, stepper);
        compute_du_avg(sim_units[su_id]->discretization, stepper);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            auto& boundary = dbound.data.boundary[dbound.bound_id];

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

        compute_dze(sim_units[su_id]->discretization, stepper);
        compute_du(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }
    /* First derivatives end */

    /* Second derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        compute_ddu_avg(sim_units[su_id]->discretization, stepper);

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

        compute_ddu(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }
    /* Second derivatives end */
}
}
}

#endif