#ifndef SWE_CS_SL_OMPI_HPP
#define SWE_CS_SL_OMPI_HPP

// Implementation of Cockburn-Shu slope limiter
#include "swe_CS_slope_limiter.hpp"
#include "swe_trouble_check.hpp"

namespace SWE {
template <typename StepperType, typename OMPISimUnitType>
void CS_slope_limiter_ompi(StepperType& stepper,
                           std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units,
                           const uint begin_sim_id,
                           const uint end_sim_id,
                           uint comm_type) {
    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Exchanging slope limiting data" << std::endl;
        }

        sim_units[su_id]->communicator.ReceiveAll(comm_type, stepper.GetTimestamp());

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { slope_limiting_prepare_element_kernel(stepper, elt); });

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper, comm_type](auto& dbound) {
            slope_limiting_distributed_boundary_send_kernel(stepper, dbound, comm_type);
        });

        sim_units[su_id]->communicator.SendAll(comm_type, stepper.GetTimestamp());
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting slope limiting work before receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachInterface(
            [&stepper](auto& intface) { slope_limiting_prepare_interface_kernel(stepper, intface); });

        sim_units[su_id]->discretization.mesh.CallForEachBoundary(
            [&stepper](auto& bound) { slope_limiting_prepare_boundary_kernel(stepper, bound); });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished slope limiting work before receive" << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile()
                << "Starting to wait on slope limiting receive with timestamp: " << stepper.GetTimestamp() << std::endl;
        }

        sim_units[su_id]->communicator.WaitAllReceives(comm_type, stepper.GetTimestamp());

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Starting slope limiting work after receive" << std::endl;
        }

        sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&stepper, comm_type](auto& dbound) {
            slope_limiting_prepare_distributed_boundary_kernel(stepper, dbound, comm_type);
        });

        check_trouble(sim_units[su_id]->discretization, stepper, comm_type);

        sim_units[su_id]->discretization.mesh.CallForEachElement(
            [&stepper](auto& elt) { slope_limiting_kernel(stepper, elt); });

        if (sim_units[su_id]->writer.WritingVerboseLog()) {
            sim_units[su_id]->writer.GetLogFile() << "Finished slope limiting work after receive" << std::endl
                                                  << std::endl;
        }
    }

    for (uint su_id = begin_sim_id; su_id < end_sim_id; ++su_id) {
        sim_units[su_id]->communicator.WaitAllSends(comm_type, stepper.GetTimestamp());
    }
}
}

#endif