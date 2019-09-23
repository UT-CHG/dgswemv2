#ifndef EHDG_GN_PRE_DBATH_OMPI_HPP
#define EHDG_GN_PRE_DBATH_OMPI_HPP

#include "ehdg_gn_pre_dbath_serial.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_bathymetry_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                                  const uint begin_sim_id,
                                                  const uint end_sim_id) {
    /* First derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            // Construct message to exterior state
            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(ngp);
            for (uint gp = 0; gp < ngp; ++gp) {
                message[gp] = boundary.aux_at_gp(SWE::Auxiliaries::bath, gp);
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_dbath_rhs(sim_units[su_id]->discretization);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(ngp);
            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex            = ngp - gp - 1;
                boundary.bath_hat_at_gp[gp] = (boundary.aux_at_gp(SWE::Auxiliaries::bath, gp) + message[gp_ex]) / 2.0;
            }

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dbath, dir) +=
                    dbound.IntegrationPhi(vec_cw_mult(boundary.bath_hat_at_gp, row(dbound.surface_normal, dir)));
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state = elt.data.state[0];
            state.dbath = elt.ApplyMinv(state.dbath);
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
    /* First derivatives end */

    /* Second derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            boundary.dbath_hat_at_gp = dbound.ComputeUgp(state.dbath);

            // Construct message to exterior state
            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_dimensions * ngp);
            for (uint gp = 0; gp < ngp; ++gp) {
                for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                    message[GN::n_dimensions * gp + dbath] = boundary.dbath_hat_at_gp(dbath, gp);
                }
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_ddbath_rhs(sim_units[su_id]->discretization);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_dimensions * ngp);
            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex = ngp - gp - 1;
                for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                    boundary.dbath_hat_at_gp(dbath, gp) =
                        (boundary.dbath_hat_at_gp(dbath, gp) + message[GN::n_dimensions * gp_ex + dbath]) / 2.0;
                }
            }

            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddbath, GN::n_dimensions * dbath + dir) += dbound.IntegrationPhi(
                        vec_cw_mult(row(boundary.dbath_hat_at_gp, dbath), row(dbound.surface_normal, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state  = elt.data.state[0];
            state.ddbath = elt.ApplyMinv(state.ddbath);
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
    /* Second derivatives end */

    /* Trird derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            boundary.ddbath_hat_at_gp = dbound.ComputeUgp(state.ddbath);

            // Construct message to exterior state
            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_ddbath_terms * ngp);
            for (uint gp = 0; gp < ngp; ++gp) {
                for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                    message[GN::n_ddbath_terms * gp + ddbath] = boundary.ddbath_hat_at_gp(ddbath, gp);
                }
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_dddbath_rhs(sim_units[su_id]->discretization);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            const uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_ddbath_terms * ngp);
            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);
            for (uint gp = 0; gp < ngp; ++gp) {
                const uint gp_ex = ngp - gp - 1;
                for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                    boundary.ddbath_hat_at_gp(ddbath, gp) =
                        (boundary.ddbath_hat_at_gp(ddbath, gp) + message[GN::n_ddbath_terms * gp_ex + ddbath]) / 2.0;
                }
            }

            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.dddbath, GN::n_dimensions * ddbath + dir) += dbound.IntegrationPhi(
                        vec_cw_mult(row(boundary.ddbath_hat_at_gp, ddbath), row(dbound.surface_normal, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state            = elt.data.state[0];
            auto& internal         = elt.data.internal;
            state.dddbath          = elt.ApplyMinv(state.dddbath);
            internal.dddbath_at_gp = elt.ComputeUgp(state.dddbath);
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
    /* Third derivatives end */
}
}
}

#endif