#ifndef EHDG_GN_PROC_DERIVATIVES_OMPI_HPP
#define EHDG_GN_PROC_DERIVATIVES_OMPI_HPP

#include "general_definitions.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
/* First derivatives begin */
#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = dbound.data.state[stage];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            boundary.q_at_gp = dbound.ComputeUgp(state.q);

            row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

            boundary.ze_hat_at_gp = row(boundary.q_at_gp, SWE::Variables::ze);

            row(boundary.u_hat_at_gp, GlobalCoord::x) =
                vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
            row(boundary.u_hat_at_gp, GlobalCoord::y) =
                vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

            // Construct message to exterior state
            std::vector<double> message;

            message.reserve(SWE::n_variables * dbound.data.get_ngp_boundary(dbound.bound_id));

            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                message.push_back(boundary.ze_hat_at_gp[gp]);
                message.push_back(boundary.u_hat_at_gp(GlobalCoord::x, gp));
                message.push_back(boundary.u_hat_at_gp(GlobalCoord::y, gp));
            }

            // Set message to send buffer
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, sim_units[su_id]->stepper.GetTimestamp());
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = elt.data.state[stage];
            auto& internal = elt.data.internal;

            internal.q_at_gp = elt.ComputeUgp(state.q);

            row(internal.aux_at_gp, SWE::Auxiliaries::h) =
                row(internal.q_at_gp, SWE::Variables::ze) + row(internal.aux_at_gp, SWE::Auxiliaries::bath);

            row(internal.u_at_gp, GlobalCoord::x) =
                vec_cw_div(row(internal.q_at_gp, SWE::Variables::qx), row(internal.aux_at_gp, SWE::Auxiliaries::h));
            row(internal.u_at_gp, GlobalCoord::y) =
                vec_cw_div(row(internal.q_at_gp, SWE::Variables::qy), row(internal.aux_at_gp, SWE::Auxiliaries::h));

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dze, dir) = -elt.IntegrationDPhi(dir, row(internal.q_at_gp, SWE::Variables::ze));
            }

            for (uint u = 0; u < GN::n_dimensions; ++u) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.du, GN::n_dimensions * u + dir) = -elt.IntegrationDPhi(dir, row(internal.u_at_gp, u));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state_in    = intface.data_in.state[stage];
            auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

            auto& state_ex    = intface.data_ex.state[stage];
            auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

            boundary_in.q_at_gp = intface.ComputeUgpIN(state_in.q);
            boundary_ex.q_at_gp = intface.ComputeUgpEX(state_ex.q);

            row(boundary_in.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary_in.q_at_gp, SWE::Variables::ze) + row(boundary_in.aux_at_gp, SWE::Auxiliaries::bath);

            row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary_ex.q_at_gp, SWE::Variables::ze) + row(boundary_ex.aux_at_gp, SWE::Auxiliaries::bath);

            row(boundary_in.u_hat_at_gp, GlobalCoord::x) = vec_cw_div(row(boundary_in.q_at_gp, SWE::Variables::qx),
                                                                      row(boundary_in.aux_at_gp, SWE::Auxiliaries::h));
            row(boundary_in.u_hat_at_gp, GlobalCoord::y) = vec_cw_div(row(boundary_in.q_at_gp, SWE::Variables::qy),
                                                                      row(boundary_in.aux_at_gp, SWE::Auxiliaries::h));

            row(boundary_ex.u_hat_at_gp, GlobalCoord::x) = vec_cw_div(row(boundary_ex.q_at_gp, SWE::Variables::qx),
                                                                      row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h));
            row(boundary_ex.u_hat_at_gp, GlobalCoord::y) = vec_cw_div(row(boundary_ex.q_at_gp, SWE::Variables::qy),
                                                                      row(boundary_ex.aux_at_gp, SWE::Auxiliaries::h));

            uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

                boundary_in.ze_hat_at_gp[gp] =
                    (boundary_in.q_at_gp(SWE::Variables::ze, gp) + boundary_ex.q_at_gp(SWE::Variables::ze, gp_ex)) /
                    2.0;

                boundary_ex.ze_hat_at_gp[gp_ex] = boundary_in.ze_hat_at_gp[gp];

                for (uint u = 0; u < GN::n_dimensions; ++u) {
                    boundary_in.u_hat_at_gp(u, gp) =
                        (boundary_in.u_hat_at_gp(u, gp) + boundary_ex.u_hat_at_gp(u, gp_ex)) / 2.0;

                    boundary_ex.u_hat_at_gp(u, gp_ex) = boundary_in.u_hat_at_gp(u, gp);
                }
            }

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.dze, dir) += intface.IntegrationPhiIN(
                    vec_cw_mult(boundary_in.ze_hat_at_gp, row(intface.surface_normal_in, dir)));

                row(state_ex.dze, dir) += intface.IntegrationPhiEX(
                    vec_cw_mult(boundary_ex.ze_hat_at_gp, row(intface.surface_normal_ex, dir)));
            }

            for (uint u = 0; u < GN::n_dimensions; ++u) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state_in.du, GN::n_dimensions * u + dir) += intface.IntegrationPhiIN(
                        vec_cw_mult(row(boundary_in.u_hat_at_gp, u), row(intface.surface_normal_in, dir)));

                    row(state_ex.du, GN::n_dimensions * u + dir) += intface.IntegrationPhiEX(
                        vec_cw_mult(row(boundary_ex.u_hat_at_gp, u), row(intface.surface_normal_ex, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = bound.data.state[stage];
            auto& boundary = bound.data.boundary[bound.bound_id];

            boundary.q_at_gp = bound.ComputeUgp(state.q);

            row(boundary.aux_at_gp, SWE::Auxiliaries::h) =
                row(boundary.q_at_gp, SWE::Variables::ze) + row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

            boundary.ze_hat_at_gp = row(boundary.q_at_gp, SWE::Variables::ze);

            row(boundary.u_hat_at_gp, GlobalCoord::x) =
                vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qx), row(boundary.aux_at_gp, SWE::Auxiliaries::h));
            row(boundary.u_hat_at_gp, GlobalCoord::y) =
                vec_cw_div(row(boundary.q_at_gp, SWE::Variables::qy), row(boundary.aux_at_gp, SWE::Auxiliaries::h));

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dze, dir) +=
                    bound.IntegrationPhi(vec_cw_mult(boundary.ze_hat_at_gp, row(bound.surface_normal, dir)));
            }

            for (uint u = 0; u < GN::n_dimensions; ++u) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.du, GN::n_dimensions * u + dir) +=
                        bound.IntegrationPhi(vec_cw_mult(row(boundary.u_hat_at_gp, u), row(bound.surface_normal, dir)));
                }
            }
        });
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = dbound.data.state[stage];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            std::vector<double> message;

            uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

            message.resize(SWE::n_variables * ngp);

            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);

            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

                boundary.ze_hat_at_gp[gp] = (boundary.ze_hat_at_gp[gp] + message[SWE::n_variables * gp_ex]) / 2.0;
                boundary.u_hat_at_gp(GlobalCoord::x, gp) =
                    (boundary.u_hat_at_gp(GlobalCoord::x, gp) + message[SWE::n_variables * gp_ex + 1]) / 2.0;
                boundary.u_hat_at_gp(GlobalCoord::y, gp) =
                    (boundary.u_hat_at_gp(GlobalCoord::y, gp) + message[SWE::n_variables * gp_ex + 2]) / 2.0;
            }

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dze, dir) +=
                    dbound.IntegrationPhi(vec_cw_mult(boundary.ze_hat_at_gp, row(dbound.surface_normal, dir)));
            }

            for (uint u = 0; u < GN::n_dimensions; ++u) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.du, GN::n_dimensions * u + dir) += dbound.IntegrationPhi(
                        vec_cw_mult(row(boundary.u_hat_at_gp, u), row(dbound.surface_normal, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state = elt.data.state[stage];

            state.dze = elt.ApplyMinv(state.dze);
            state.du  = elt.ApplyMinv(state.du);
        });
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, sim_units[su_id]->stepper.GetTimestamp());
    }
/* First derivatives end */

/* Second derivatives begin */
#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::derivatives, sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = dbound.data.state[stage];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            boundary.du_hat_at_gp = dbound.ComputeUgp(state.du);

            // Construct message to exterior state
            std::vector<double> message;

            message.reserve(GN::n_du_terms * dbound.data.get_ngp_boundary(dbound.bound_id));

            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                for (uint du = 0; du < GN::n_du_terms; ++du) {
                    message.push_back(boundary.du_hat_at_gp(du, gp));
                }
            }

            // Set message to send buffer
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::derivatives, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::derivatives, sim_units[su_id]->stepper.GetTimestamp());
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = elt.data.state[stage];
            auto& internal = elt.data.internal;

            internal.du_at_gp = elt.ComputeUgp(state.du);

            for (uint du = 0; du < GN::n_du_terms; ++du) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddu, GN::n_dimensions * du + dir) = -elt.IntegrationDPhi(dir, row(internal.du_at_gp, du));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

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
                        vec_cw_mult(row(boundary_in.du_hat_at_gp, du), row(intface.surface_normal_in, dir)));

                    row(state_ex.ddu, GN::n_dimensions * du + dir) += intface.IntegrationPhiEX(
                        vec_cw_mult(row(boundary_ex.du_hat_at_gp, du), row(intface.surface_normal_ex, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = bound.data.state[stage];
            auto& boundary = bound.data.boundary[bound.bound_id];

            boundary.du_hat_at_gp = bound.ComputeUgp(state.du);

            for (uint du = 0; du < GN::n_du_terms; ++du) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddu, GN::n_dimensions * du + dir) += bound.IntegrationPhi(
                        vec_cw_mult(row(boundary.du_hat_at_gp, du), row(bound.surface_normal, dir)));
                }
            }
        });
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::derivatives,
                                                       sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = dbound.data.state[stage];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            std::vector<double> message;

            uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

            message.resize(GN::n_du_terms * ngp);

            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::derivatives, message);

            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

                for (uint du = 0; du < GN::n_du_terms; ++du) {
                    boundary.du_hat_at_gp(du, gp) =
                        (boundary.du_hat_at_gp(du, gp) + message[GN::n_du_terms * gp_ex + du]) / 2.0;
                }
            }

            for (uint du = 0; du < GN::n_du_terms; ++du) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddu, GN::n_dimensions * du + dir) += dbound.IntegrationPhi(
                        vec_cw_mult(row(boundary.du_hat_at_gp, du), row(dbound.surface_normal, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state    = elt.data.state[stage];
            auto& internal = elt.data.internal;

            state.ddu = elt.ApplyMinv(state.ddu);

            internal.ddu_at_gp = elt.ComputeUgp(state.ddu);
        });
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::derivatives, sim_units[su_id]->stepper.GetTimestamp());
    }
    /* Second derivatives end */
}
}
}

#endif