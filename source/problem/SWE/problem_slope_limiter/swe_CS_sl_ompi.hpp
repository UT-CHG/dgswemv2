#ifndef SWE_CS_SL_OMPI_HPP
#define SWE_CS_SL_OMPI_HPP

// Implementation of Cockburn-Shu slope limiter
#include "swe_CS_slope_limiter.hpp"

namespace SWE {
template <typename OMPISimUnitType>
void CS_slope_limiter_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                           uint begin_sim_id,
                           uint end_sim_id,
                           uint comm_type) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Exchanging slope limiting data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(comm_type, sim_units[su_id]->stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { slope_limiting_prepare_element_kernel(sim_units[su_id]->stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [&sim_units, su_id, comm_type](auto& dbound) {
                slope_limiting_distributed_boundary_send_kernel(sim_units[su_id]->stepper, dbound, comm_type);
            });

        sim_units[su_id]->communicator.SendAll(comm_type, sim_units[su_id]->stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting slope limiting work before receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
            slope_limiting_prepare_interface_kernel(sim_units[su_id]->stepper, intface);
        });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary([&sim_units, su_id](auto& bound) {
            slope_limiting_prepare_boundary_kernel(sim_units[su_id]->stepper, bound);
        });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished slope limiting work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting to wait on slope limiting receive with timestamp: "
                                                  << sim_units[su_id]->stepper.GetTimestamp() << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(comm_type, sim_units[su_id]->stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting slope limiting work after receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary(
            [&sim_units, su_id, comm_type](auto& dbound) {
                slope_limiting_prepare_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound, comm_type);
            });

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&sim_units, su_id](auto& elt) { slope_limiting_kernel(sim_units[su_id]->stepper, elt); });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished slope limiting work after receive" << std::endl
                                                  << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; su_id++) {
        sim_units[su_id]->communicator.WaitAllSends(comm_type, sim_units[su_id]->stepper.GetTimestamp());
    }
}
}

#endif