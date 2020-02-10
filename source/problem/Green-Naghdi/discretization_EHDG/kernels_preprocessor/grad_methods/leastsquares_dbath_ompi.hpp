#ifndef LEASTSQUARES_DBATH_OMPI_HPP
#define LEASTSQUARES_DBATH_OMPI_HPP

#include "leastsquares_dbath_serial.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_bathymetry_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                                  ProblemGlobalDataType& global_data,
                                                  const uint begin_sim_id,
                                                  const uint end_sim_id,
                                                  const uint stage) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachElement([stage](auto& elt) {
            auto& state                = elt.data.state[stage];
            auto& derivative           = elt.data.derivative;
            derivative.bath_lin        = elt.ProjectBasisToLinear(row(state.aux, SWE::Auxiliaries::bath));
            derivative.bath_at_baryctr = elt.ComputeLinearUbaryctr(derivative.bath_lin);
            derivative.bath_at_midpts  = elt.ComputeLinearUmidpts(derivative.bath_lin);
        });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(1, derivative.bath_at_baryctr);
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        compute_dbath_ls(sim_units[su_id]->discretization, stage);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }

#if defined(B_RECON_INT) || defined(B_RECON_LS)
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(GN::n_dimensions);
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                message[dbath] = derivative.dbath_at_baryctr[dbath];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        reconstruct_dbath(sim_units[su_id]->discretization, global_data, stage);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
#elif defined(B_RECON_AVG)
    reconstruct_dbath(sim_units, global_data, stage);
#endif

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& derivative            = elt.data.derivative;
            derivative.dbath_at_baryctr = elt.ComputeLinearUbaryctr(derivative.dbath_lin);
            derivative.dbath_at_midpts  = elt.ComputeLinearUmidpts(derivative.dbath_lin);
        });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(GN::n_dimensions);
            for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                message[dbath] = derivative.dbath_at_baryctr[dbath];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        compute_ddbath_ls(sim_units[su_id]->discretization, stage);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }

#if defined(B_RECON_INT) || defined(B_RECON_LS)
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(GN::n_ddbath_terms);
            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                message[ddbath] = derivative.ddbath_at_baryctr[ddbath];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        reconstruct_ddbath(sim_units[su_id]->discretization, global_data, stage);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
#elif defined(B_RECON_AVG)
    reconstruct_ddbath(sim_units, global_data, stage);
#endif
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachElement([](auto& elt) {
            auto& derivative             = elt.data.derivative;
            derivative.ddbath_at_baryctr = elt.ComputeLinearUbaryctr(derivative.ddbath_lin);
            derivative.ddbath_at_midpts  = elt.ComputeLinearUmidpts(derivative.ddbath_lin);
        });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            auto& boundary   = dbound.data.boundary[dbound.bound_id];
            const uint ngp   = dbound.data.get_ngp_boundary(dbound.bound_id);
            std::vector<double> message(GN::n_ddbath_terms + GN::n_dimensions * ngp);
            for (uint ddbath = 0; ddbath < GN::n_ddbath_terms; ++ddbath) {
                message[ddbath] = derivative.ddbath_at_baryctr[ddbath];
            }
            for (uint gp = 0; gp < ngp; ++gp) {
                for (uint dbath = 0; dbath < GN::n_dimensions; ++dbath) {
                    message[GN::n_ddbath_terms + GN::n_dimensions * gp + dbath] = boundary.dbath_hat_at_gp(dbath, gp);
                }
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        compute_dddbath_ls(sim_units[su_id]->discretization, stage);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }

#if defined(B_RECON_INT) || defined(B_RECON_LS)
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& derivative = dbound.data.derivative;
            std::vector<double> message(GN::n_dddbath_terms);
            for (uint dddbath = 0; dddbath < GN::n_dddbath_terms; ++dddbath) {
                message[dddbath] = derivative.dddbath_at_baryctr[dddbath];
            }
            dbound.boundary_condition.exchanger.SetToSendBuffer(CommTypes::dbath, message);
        });

        sim_units[su_id]->communicator.SendAll(CommTypes::dbath, 0);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        reconstruct_dddbath(sim_units[su_id]->discretization, global_data, stage);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
#elif defined(B_RECON_AVG)
    reconstruct_dddbath(sim_units, global_data, stage);
#endif
}
}
}

#endif