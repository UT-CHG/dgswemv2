#ifndef IHDG_SWE_PROC_INIT_ITER_HPP
#define IHDG_SWE_PROC_INIT_ITER_HPP

namespace SWE {
namespace IHDG {
template <typename StepperType, typename ProblemType>
void Problem::init_iteration(const StepperType& stepper, HDGDiscretization<ProblemType>& discretization) {
    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::init_volume_kernel(stepper, elt); });

    discretization.mesh.CallForEachElement([&stepper](auto& elt) { Problem::init_source_kernel(stepper, elt); });

    discretization.mesh.CallForEachInterface(
        [&stepper](auto& intface) { Problem::init_interface_kernel(stepper, intface); });

    discretization.mesh.CallForEachBoundary([&stepper](auto& bound) { Problem::init_boundary_kernel(stepper, bound); });

    discretization.mesh.CallForEachDistributedBoundary(
        [&stepper](auto& bound) { Problem::init_distributed_boundary_kernel(stepper, bound); });

    discretization.mesh_skeleton.CallForEachEdgeInterface(
        [&stepper](auto& edge_int) { Problem::init_edge_interface_kernel(stepper, edge_int); });

    discretization.mesh_skeleton.CallForEachEdgeBoundary(
        [&stepper](auto& edge_bound) { Problem::init_edge_boundary_kernel(stepper, edge_bound); });

    discretization.mesh_skeleton.CallForEachEdgeDistributed(
        [&stepper](auto& edge_bound) { Problem::init_edge_distributed_kernel(stepper, edge_bound); });
}
}
}

#endif