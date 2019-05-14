#ifndef EHDG_SWE_PROC_HPX_STAGE_HPP
#define EHDG_SWE_PROC_HPX_STAGE_HPP

#include "ehdg_swe_kernels_processor.hpp"

namespace SWE {
namespace EHDG {
template <typename HPXSimUnitType>
auto Problem::stage_hpx(HPXSimUnitType* sim_unit) {
    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Current (time, stage): (" << sim_unit->stepper.GetTimeAtCurrentStage() << ','
                                      << sim_unit->stepper.GetStage() << ')' << std::endl;

        sim_unit->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    hpx::future<void> receive_future =
        sim_unit->communicator.ReceiveAll(CommTypes::bound_state, sim_unit->stepper.GetTimestamp());

    sim_unit->discretization.mesh.CallForEachDistributedBoundary(
        [sim_unit](auto& dbound) { Problem::global_distributed_boundary_kernel(sim_unit->stepper, dbound); });

    sim_unit->communicator.SendAll(CommTypes::bound_state, sim_unit->stepper.GetTimestamp());

    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    /* Global Pre Receive Step */
    sim_unit->discretization.mesh.CallForEachInterface(
        [sim_unit](auto& intface) { Problem::global_interface_kernel(sim_unit->stepper, intface); });

    sim_unit->discretization.mesh.CallForEachBoundary(
        [sim_unit](auto& bound) { Problem::global_boundary_kernel(sim_unit->stepper, bound); });

    sim_unit->discretization.mesh_skeleton.CallForEachEdgeInterface(
        [sim_unit](auto& edge_int) { edge_int.interface.specialization.ComputeNumericalFlux(edge_int); });

    sim_unit->discretization.mesh_skeleton.CallForEachEdgeBoundary([sim_unit](auto& edge_bound) {
        edge_bound.boundary.boundary_condition.ComputeNumericalFlux(sim_unit->stepper, edge_bound);
    });

    /* Global Pre Receive Step */

    /* Local Pre Receive Step */
    sim_unit->discretization.mesh.CallForEachElement(
        [sim_unit](auto& elt) { Problem::local_volume_kernel(sim_unit->stepper, elt); });

    sim_unit->discretization.mesh.CallForEachElement(
        [sim_unit](auto& elt) { Problem::local_source_kernel(sim_unit->stepper, elt); });

    sim_unit->discretization.mesh.CallForEachInterface(
        [sim_unit](auto& intface) { Problem::local_interface_kernel(sim_unit->stepper, intface); });

    sim_unit->discretization.mesh.CallForEachBoundary(
        [sim_unit](auto& bound) { Problem::local_boundary_kernel(sim_unit->stepper, bound); });
    /* Local Pre Receive Step */

    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Finished work before receive" << std::endl
                                      << "Starting to wait on receive with timestamp: "
                                      << sim_unit->stepper.GetTimestamp() << std::endl;
    }

    hpx::future<void> stage_future = receive_future.then([sim_unit](auto&&) {
        if (sim_unit->writer.WritingVerboseLog()) {
            sim_unit->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        /* Global Post Receive Step */
        sim_unit->discretization.mesh_skeleton.CallForEachEdgeDistributed([sim_unit](auto& edge_dbound) {
            edge_dbound.boundary.boundary_condition.ComputeNumericalFlux(edge_dbound);
        });
        /* Global Post Receive Step */

        /* Local Post Receive Step */
        sim_unit->discretization.mesh.CallForEachDistributedBoundary(
            [sim_unit](auto& dbound) { Problem::local_distributed_boundary_kernel(sim_unit->stepper, dbound); });

        UpdateKernel update_kernel(sim_unit->stepper);
        sim_unit->discretization.mesh.CallForEachElement(update_kernel);

        ++(sim_unit->stepper);
        /* Local Post Receive Step */

        if (sim_unit->writer.WritingVerboseLog()) {
            sim_unit->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    });

    if (SWE::PostProcessing::wetting_drying) {
        stage_future = stage_future.then([sim_unit](auto&& f) {
            f.get();  // check for exceptions

            sim_unit->discretization.mesh.CallForEachElement(
                [sim_unit](auto& elt) { wetting_drying_kernel(sim_unit->stepper, elt); });
        });
    }

    if (SWE::PostProcessing::slope_limiting) {
        stage_future = stage_future.then([sim_unit](auto&& f) {
            f.get();  // check for exceptions

            return CS_slope_limiter_hpx(sim_unit, CommTypes::baryctr_state);
        });
    }

    return stage_future.then([sim_unit](auto&& f) {
        f.get();  // check for exceptions

        sim_unit->discretization.mesh.CallForEachElement([sim_unit](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(sim_unit->stepper, elt);

            if (nan_found)
                hpx::terminate();
        });
    });
}
}
}

#endif