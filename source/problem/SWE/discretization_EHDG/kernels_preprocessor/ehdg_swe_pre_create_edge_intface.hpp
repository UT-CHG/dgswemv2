#ifndef EHDG_SWE_PRE_CREATE_EDGE_INTFACE_HPP
#define EHDG_SWE_PRE_CREATE_EDGE_INTFACE_HPP

namespace SWE {
namespace EHDG {
void Problem::create_edge_interfaces(ProblemMeshType& mesh,
                                     ProblemMeshSkeletonType& mesh_skeleton,
                                     Writer<Problem>& writer) {
    using InterfaceTypes = Geometry::InterfaceTypeTuple<Data, IS::Internal>;

    using EdgeInterfaceTypes =
        Geometry::EdgeInterfaceTypeTuple<EdgeData, Geometry::InterfaceTypeTuple<Data, IS::Internal>>::Type;

    uint old_edge_internal = mesh_skeleton.GetNumberEdgeInterfaces();

    using InterfaceTypeInternal = typename std::tuple_element<0, InterfaceTypes>::type;
    mesh.CallForEachInterfaceOfType<InterfaceTypeInternal>([&mesh_skeleton, &writer](auto& intface) {
        using EdgeInterfaceTypeInternal = typename std::tuple_element<0, EdgeInterfaceTypes>::type;

        mesh_skeleton.template CreateEdgeInterface<EdgeInterfaceTypeInternal>(intface);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of internal edges: "
                            << mesh_skeleton.GetNumberEdgeInterfaces() - old_edge_internal << std::endl;
    }
}
}
}

#endif
