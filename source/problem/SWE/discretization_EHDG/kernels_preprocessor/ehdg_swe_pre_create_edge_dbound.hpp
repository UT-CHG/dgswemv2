#ifndef EHDG_SWE_PRE_CREATE_EDGE_DBOUND_HPP
#define EHDG_SWE_PRE_CREATE_EDGE_DBOUND_HPP

namespace SWE {
namespace EHDG {
void Problem::create_edge_distributeds(ProblemMeshType& mesh,
                                       ProblemMeshSkeletonType& mesh_skeleton,
                                       ProblemWriterType& writer) {
    using DistributedBoundaryTypes = Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>;

    using EdgeDistributedTypes =
        Geometry::EdgeDistributedTypeTuple<EdgeData,
                                           Geometry::DistributedBoundaryTypeTuple<Data, DBC::Distributed>>::Type;

    uint old_edge_distributed = mesh_skeleton.GetNumberEdgeDistributeds();

    using DBTypeDistributed = typename std::tuple_element<0, DistributedBoundaryTypes>::type;
    mesh.CallForEachDistributedBoundaryOfType<DBTypeDistributed>([&mesh_skeleton, &writer](auto& dbound) {
        using EDTypeDistributed = typename std::tuple_element<0, EdgeDistributedTypes>::type;

        uint locality_in = dbound.boundary_condition.exchanger.locality_in;
        uint submesh_in  = dbound.boundary_condition.exchanger.submesh_in;

        uint locality_ex = dbound.boundary_condition.exchanger.locality_ex;
        uint submesh_ex  = dbound.boundary_condition.exchanger.submesh_ex;

        bool ccw;

        if (locality_in < locality_ex || (locality_in == locality_ex && submesh_in < submesh_ex)) {
            ccw = true;  // trace basis is defined counterclockwise around element
        } else if (locality_in > locality_ex || (locality_in == locality_ex && submesh_in > submesh_ex)) {
            ccw = false;  // trace basis is defined clockwise around element
        }

        mesh_skeleton.template CreateEdgeDistributed<EDTypeDistributed>(dbound, ccw);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of distributed edges: "
                            << mesh_skeleton.GetNumberEdgeDistributeds() - old_edge_distributed << std::endl;
    }
}
}
}

#endif
