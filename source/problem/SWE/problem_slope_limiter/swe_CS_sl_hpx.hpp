#ifndef SWE_CS_SL_HPX_HPP
#define SWE_CS_SL_HPX_HPP

// Implementation of Cockburn-Shu slope limiter
#include "swe_CS_slope_limiter.hpp"
#include "swe_trouble_check.hpp"

namespace SWE {
template <typename HPXSimUnitType>
auto CS_slope_limiter_hpx(HPXSimUnitType* sim_unit, uint comm_type) {
    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Exchanging slope limiting data" << std::endl;
    }

    hpx::future<void> receive_future = sim_unit->communicator.ReceiveAll(comm_type, sim_unit->stepper.GetTimestamp());

    sim_unit->discretization.mesh.CallForEachElement(
        [sim_unit](auto& elt) { slope_limiting_prepare_element_kernel(sim_unit->stepper, elt); });

    sim_unit->discretization.mesh.CallForEachDistributedBoundary([sim_unit, comm_type](auto& dbound) {
        slope_limiting_distributed_boundary_send_kernel(sim_unit->stepper, dbound, comm_type);
    });

    sim_unit->communicator.SendAll(comm_type, sim_unit->stepper.GetTimestamp());

    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Starting slope limiting work before receive" << std::endl;
    }

    sim_unit->discretization.mesh.CallForEachInterface(
        [sim_unit](auto& intface) { slope_limiting_prepare_interface_kernel(sim_unit->stepper, intface); });

    sim_unit->discretization.mesh.CallForEachBoundary(
        [sim_unit](auto& bound) { slope_limiting_prepare_boundary_kernel(sim_unit->stepper, bound); });

    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Finished slope limiting work before receive" << std::endl
                                      << "Starting to wait on slope limiting receive with timestamp: "
                                      << sim_unit->stepper.GetTimestamp() << std::endl;
    }

    return receive_future.then([sim_unit, comm_type](auto&&) {
        if (sim_unit->writer.WritingVerboseLog()) {
            sim_unit->writer.GetLogFile() << "Starting slope limiting work after receive" << std::endl;
        }

        sim_unit->discretization.mesh.CallForEachDistributedBoundary([sim_unit, comm_type](auto& dbound) {
            slope_limiting_prepare_distributed_boundary_kernel(sim_unit->stepper, dbound, comm_type);
        });

        check_trouble(sim_unit->discretization, sim_unit->stepper, comm_type);

        sim_unit->discretization.mesh.CallForEachElement(
            [sim_unit](auto& elt) { slope_limiting_kernel(sim_unit->stepper, elt); });

        if (sim_unit->writer.WritingVerboseLog()) {
            sim_unit->writer.GetLogFile() << "Finished slope limiting work after receive" << std::endl << std::endl;
        }
    });
}
}

#endif