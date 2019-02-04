#ifndef EHDG_GN_PRE_DBATH_OMPI_HPP
#define EHDG_GN_PRE_DBATH_OMPI_HPP

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_bathymetry_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                                  uint begin_sim_id,
                                                  uint end_sim_id) {
    /* First derivatives begin */
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            // Construct message to exterior state
            std::vector<double> message;

            message.reserve(dbound.data.get_ngp_boundary(dbound.bound_id));

            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                message.push_back(boundary.aux_at_gp(SWE::Auxiliaries::bath, gp));
            }

            // Set message to send buffer
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state    = elt.data.state[0];
            auto& internal = elt.data.internal;

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dbath, dir) = -elt.IntegrationDPhi(dir, row(internal.aux_at_gp, SWE::Auxiliaries::bath));
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([](auto& intface) {
            auto& state_in    = intface.data_in.state[0];
            auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

            auto& state_ex    = intface.data_ex.state[0];
            auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

            uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

                boundary_in.bath_hat_at_gp[gp] = (boundary_in.aux_at_gp(SWE::Auxiliaries::bath, gp) +
                                                  boundary_ex.aux_at_gp(SWE::Auxiliaries::bath, gp_ex)) /
                                                 2.0;

                boundary_ex.bath_hat_at_gp[gp_ex] = boundary_in.bath_hat_at_gp[gp];
            }

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state_in.dbath, dir) += intface.IntegrationPhiIN(
                    vec_cw_mult(boundary_in.bath_hat_at_gp, row(intface.surface_normal_in, dir)));

                row(state_ex.dbath, dir) += intface.IntegrationPhiEX(
                    vec_cw_mult(boundary_ex.bath_hat_at_gp, row(intface.surface_normal_ex, dir)));
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([](auto& bound) {
            auto& state    = bound.data.state[0];
            auto& boundary = bound.data.boundary[bound.bound_id];

            boundary.bath_hat_at_gp = row(boundary.aux_at_gp, SWE::Auxiliaries::bath);

            for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                row(state.dbath, dir) +=
                    bound.IntegrationPhi(vec_cw_mult(boundary.bath_hat_at_gp, row(bound.surface_normal, dir)));
            }
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            std::vector<double> message;

            uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

            message.resize(ngp);

            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);

            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

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
            std::vector<double> message;

            message.reserve(GN::n_dimensions * dbound.data.get_ngp_boundary(dbound.bound_id));

            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                    message.push_back(boundary.dbath_hat_at_gp(dbath, gp));
                }
            }

            // Set message to send buffer
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state    = elt.data.state[0];
            auto& internal = elt.data.internal;

            internal.dbath_at_gp = elt.ComputeUgp(state.dbath);

            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddbath, GN::n_dimensions * dbath + dir) =
                        -elt.IntegrationDPhi(dir, row(internal.dbath_at_gp, dbath));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([](auto& intface) {
            auto& state_in    = intface.data_in.state[0];
            auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

            auto& state_ex    = intface.data_ex.state[0];
            auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

            boundary_in.dbath_hat_at_gp = intface.ComputeUgpIN(state_in.dbath);
            boundary_ex.dbath_hat_at_gp = intface.ComputeUgpEX(state_ex.dbath);

            uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

                for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                    boundary_in.dbath_hat_at_gp(dbath, gp) =
                        (boundary_in.dbath_hat_at_gp(dbath, gp) + boundary_ex.dbath_hat_at_gp(dbath, gp_ex)) / 2.0;

                    boundary_ex.dbath_hat_at_gp(dbath, gp_ex) = boundary_in.dbath_hat_at_gp(dbath, gp);
                }
            }

            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state_in.ddbath, GN::n_dimensions * dbath + dir) += intface.IntegrationPhiIN(
                        vec_cw_mult(row(boundary_in.dbath_hat_at_gp, dbath), row(intface.surface_normal_in, dir)));

                    row(state_ex.ddbath, GN::n_dimensions * dbath + dir) += intface.IntegrationPhiEX(
                        vec_cw_mult(row(boundary_ex.dbath_hat_at_gp, dbath), row(intface.surface_normal_ex, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([](auto& bound) {
            auto& state    = bound.data.state[0];
            auto& boundary = bound.data.boundary[bound.bound_id];

            boundary.dbath_hat_at_gp = bound.ComputeUgp(state.dbath);

            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.ddbath, GN::n_dimensions * dbath + dir) += bound.IntegrationPhi(
                        vec_cw_mult(row(boundary.dbath_hat_at_gp, dbath), row(bound.surface_normal, dir)));
                }
            }
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            std::vector<double> message;

            uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

            message.resize(GN::n_dimensions * ngp);

            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);

            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

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
            auto& state = elt.data.state[0];

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
            std::vector<double> message;

            message.reserve(GN::n_ddbath_terms * dbound.data.get_ngp_boundary(dbound.bound_id));

            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                    message.push_back(boundary.ddbath_hat_at_gp(ddbath, gp));
                }
            }

            // Set message to send buffer
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& state    = elt.data.state[0];
            auto& internal = elt.data.internal;

            internal.ddbath_at_gp = elt.ComputeUgp(state.ddbath);

            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.dddbath, GN::n_dimensions * ddbath + dir) =
                        -elt.IntegrationDPhi(dir, row(internal.ddbath_at_gp, ddbath));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachInterface([](auto& intface) {
            auto& state_in    = intface.data_in.state[0];
            auto& boundary_in = intface.data_in.boundary[intface.bound_id_in];

            auto& state_ex    = intface.data_ex.state[0];
            auto& boundary_ex = intface.data_ex.boundary[intface.bound_id_ex];

            boundary_in.ddbath_hat_at_gp = intface.ComputeUgpIN(state_in.ddbath);
            boundary_ex.ddbath_hat_at_gp = intface.ComputeUgpEX(state_ex.ddbath);

            uint ngp = intface.data_in.get_ngp_boundary(intface.bound_id_in);
            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

                for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                    boundary_in.ddbath_hat_at_gp(ddbath, gp) =
                        (boundary_in.ddbath_hat_at_gp(ddbath, gp) + boundary_ex.ddbath_hat_at_gp(ddbath, gp_ex)) / 2.0;

                    boundary_ex.ddbath_hat_at_gp(ddbath, gp_ex) = boundary_in.ddbath_hat_at_gp(ddbath, gp);
                }
            }

            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state_in.dddbath, GN::n_dimensions * ddbath + dir) += intface.IntegrationPhiIN(
                        vec_cw_mult(row(boundary_in.ddbath_hat_at_gp, ddbath), row(intface.surface_normal_in, dir)));

                    row(state_ex.dddbath, GN::n_dimensions * ddbath + dir) += intface.IntegrationPhiEX(
                        vec_cw_mult(row(boundary_ex.ddbath_hat_at_gp, ddbath), row(intface.surface_normal_ex, dir)));
                }
            }
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([](auto& bound) {
            auto& state    = bound.data.state[0];
            auto& boundary = bound.data.boundary[bound.bound_id];

            boundary.ddbath_hat_at_gp = bound.ComputeUgp(state.ddbath);

            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                for (uint dir = 0; dir < GN::n_dimensions; ++dir) {
                    row(state.dddbath, GN::n_dimensions * ddbath + dir) += bound.IntegrationPhi(
                        vec_cw_mult(row(boundary.ddbath_hat_at_gp, ddbath), row(bound.surface_normal, dir)));
                }
            }
        });
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& state    = dbound.data.state[0];
            auto& boundary = dbound.data.boundary[dbound.bound_id];

            std::vector<double> message;

            uint ngp = dbound.data.get_ngp_boundary(dbound.bound_id);

            message.resize(GN::n_ddbath_terms * ngp);

            dbound.boundary_condition.exchanger.GetFromReceiveBuffer(CommTypes::dbath, message);

            uint gp_ex;
            for (uint gp = 0; gp < ngp; ++gp) {
                gp_ex = ngp - gp - 1;

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
            auto& state    = elt.data.state[0];
            auto& internal = elt.data.internal;

            state.dddbath = elt.ApplyMinv(state.dddbath);

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