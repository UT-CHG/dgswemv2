#ifndef GREENGAUSS_DBATH_OMPI_HPP
#define GREENGAUSS_DBATH_OMPI_HPP

#include "greengauss_dbath_serial.hpp"

namespace GN {
namespace EHDG {
template <typename OMPISimUnitType>
void Problem::compute_bathymetry_derivatives_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                                                  ProblemGlobalDataType& global_data,
                                                  const uint begin_sim_id,
                                                  const uint end_sim_id) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_dbath_gg(sim_units[su_id]->discretization);
    }

#if defined(D_RECON_INT) || defined(D_RECON_LS)
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

        reconstruct_dbath(sim_units[su_id]->discretization, global_data);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
#elif defined(D_RECON_AVG)
    reconstruct_dbath(sim_units, global_data);
#endif

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        compute_ddbath_gg(sim_units[su_id]->discretization);
    }

#if defined(D_RECON_INT) || defined(D_RECON_LS)
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

        reconstruct_ddbath(sim_units[su_id]->discretization, global_data);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
#elif defined(D_RECON_AVG)
    reconstruct_ddbath(sim_units, global_data);
#endif

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.ReceiveAll(CommTypes::dbath, 0);

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([](auto& dbound) {
            auto& boundary = dbound.data.boundary[dbound.bound_id];
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
        sim_units[su_id]->communicator.WaitAllReceives(CommTypes::dbath, 0);

        compute_dddbath_gg(sim_units[su_id]->discretization);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }

#if defined(D_RECON_INT) || defined(D_RECON_LS)
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

        reconstruct_dddbath(sim_units[su_id]->discretization, global_data);
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(CommTypes::dbath, 0);
    }
#elif defined(D_RECON_AVG)
    reconstruct_dddbath(sim_units, global_data);
#endif
}
}
}

#endif