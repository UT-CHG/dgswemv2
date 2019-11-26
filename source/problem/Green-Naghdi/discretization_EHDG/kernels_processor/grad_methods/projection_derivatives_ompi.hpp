#ifndef PROJECTION_DERIVATIVES_OMPI_HPP
#define PROJECTION_DERIVATIVES_OMPI_HPP

#include "projection_derivatives_serial.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                       ProblemGlobalDataType& global_data,
                                       const ESSPRKStepper& stepper,
                                       const uint begin_sim_id,
                                       const uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            const uint stage = stepper.GetStage();
            auto& state      = dbound.data.state[stage];
            auto& derivative = dbound.data.derivative;
            auto& boundary   = dbound.data.boundary[dbound.bound_id];

            boundary.q_at_gp = dbound.ComputeUgp(state.q);
            row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);
            derivative.ze_hat_at_gp[dbound.bound_id] = row(boundary.q_at_gp, SWE::Variables::ze);
            row(derivative.u_hat_at_gp[dbound.bound_id], GlobalCoord::x) =
                vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
            row(derivative.u_hat_at_gp[dbound.bound_id], GlobalCoord::y) =
                vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(SWE::n_variables * ngp);
            for (uint gp = 0; gp < ngp; ++gp) {
                message[SWE::n_variables * gp]     = derivative.ze_hat_at_gp[dbound.bound_id][gp];
                message[SWE::n_variables * gp + 1] = derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::x, gp);
                message[SWE::n_variables * gp + 2] = derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::y, gp);
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_dze_rhs(sim_units[su_id]->discretization, stepper);
        compute_du_rhs(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            const uint stage = stepper.GetStage();
            auto& state      = dbound.data.state[stage];
            auto& derivative = dbound.data.derivative;
            auto& boundary   = dbound.data.boundary[dbound.bound_id];

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(SWE::n_variables * ngp);
            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex = ngp - gp - 1;
                derivative.ze_hat_at_gp[dbound.bound_id][gp] =
                    (derivative.ze_hat_at_gp[dbound.bound_id][gp] + message[SWE::n_variables * gp_ex]) / 2.0;
                derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::x, gp) =
                    (derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::x, gp) +
                     message[SWE::n_variables * gp_ex + 1]) /
                    2.0;
                derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::y, gp) =
                    (derivative.u_hat_at_gp[dbound.bound_id](GlobalCoord::y, gp) +
                     message[SWE::n_variables * gp_ex + 2]) /
                    2.0;

                // HERE I ASSUME THAT BATHYMETRY IS CONTINUOUS
                boundary.aux_at_gp(SWE::Auxiliaries::h, gp) =
                    derivative.ze_hat_at_gp[dbound.bound_id][gp] + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
            }

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dze, dir) += dbound.IntegrationPhi(
                    vec_cw_mult(derivative.ze_hat_at_gp[dbound.bound_id], row(dbound.surface_normal, dir)));
            }

            for (uint u = 0; u < GN::n_dimensions; ++u) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.du, GN::n_dimensions * u + dir) += dbound.IntegrationPhi(
                        vec_cw_mult(row(derivative.u_hat_at_gp[dbound.bound_id], u), row(dbound.surface_normal, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            const uint stage = stepper.GetStage();
            auto& state      = elt.data.state[stage];

            state.dze = elt.ApplyMinv(state.dze);
            state.du  = elt.ApplyMinv(state.du);
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }
    /* First derivatives end */

    /* Second derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            const uint stage = stepper.GetStage();
            auto& state      = dbound.data.state[stage];
            auto& derivative = dbound.data.derivative;

            derivative.du_hat_at_gp[dbound.bound_id] = dbound.ComputeUgp(state.du);

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_du_terms * ngp);
            for (uint gp = 0; gp < ngp; ++gp) {
                for (uint du = 0; du < GN::n_du_terms; ++du) {
                    message[GN::n_du_terms * gp + du] = derivative.du_hat_at_gp[dbound.bound_id](du, gp);
                }
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_ddu_rhs(sim_units[su_id]->discretization, stepper);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper](auto& dbound) {
            const uint stage = stepper.GetStage();
            auto& state      = dbound.data.state[stage];
            auto& derivative = dbound.data.derivative;

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_du_terms * ngp);
            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex = ngp - gp - 1;
                for (uint du = 0; du < GN::n_du_terms; ++du) {
                    derivative.du_hat_at_gp[dbound.bound_id](du, gp) =
                        (derivative.du_hat_at_gp[dbound.bound_id](du, gp) + message[GN::n_du_terms * gp_ex + du]) / 2.0;
                }
            }

            for (uint du = 0; du < GN::n_du_terms; ++du) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddu, GN::n_dimensions * du + dir) += dbound.IntegrationPhi(vec_cw_mult(
                        row(derivative.du_hat_at_gp[dbound.bound_id], du), row(dbound.surface_normal, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&stepper](auto& elt) {
            const uint stage = stepper.GetStage();
            auto& state      = elt.data.state[stage];
            auto& internal   = elt.data.internal;

            state.ddu          = elt.ApplyMinv(state.ddu);
            internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, stepper.GetTimestamp());
    }
    /* Second derivatives end */
}
}
}

#endif