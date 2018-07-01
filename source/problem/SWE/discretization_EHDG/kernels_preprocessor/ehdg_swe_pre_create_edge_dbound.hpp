#ifndef EHDG_SWE_PRE_CREATE_EDGE_DBOUND_HPP
#define EHDG_SWE_PRE_CREATE_EDGE_DBOUND_HPP

namespace SWE {
namespace EHDG {
void Problem::create_edge_distributeds_kernel(ProblemMeshType& mesh,
                                              ProblemMeshSkeletonType& mesh_skeleton,
                                              Writer<Problem>& writer) {
    using DistributedBoundaryTypes = Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>;

    using EdgeDistributedTypes =
        Geometry::EdgeDistributedTypeTuple<EdgeData,
                                           Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>>::Type;

    uint old_edge_distributed = mesh_skeleton.GetNumberEdgeDistributeds();

    using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;
    mesh.CallForEachDistributedBoundaryOfType<DBTypeDistributed>([&mesh_skeleton, &writer](auto& dbound) {
        using EDTypeDistributed = typename std::tuple_element<0, EdgeDistributedTypes>::type;

        mesh_skeleton.template CreateEdgeDistributed<EDTypeDistributed>(dbound);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed edges: "
                            << mesh_skeleton.GetNumberEdgeDistributeds() - old_edge_distributed << std::endl;
    }
}
}
}

#endif
