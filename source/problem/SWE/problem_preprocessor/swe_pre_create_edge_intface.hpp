#ifndef EHDG_SWE_PRE_CREATE_EDGE_INTFACE_HPP
#define EHDG_SWE_PRE_CREATE_EDGE_INTFACE_HPP

namespace SWE {
template <typename ProblemType>
void create_edge_interfaces(typename ProblemType::ProblemMeshType& mesh,
                            typename ProblemType::ProblemMeshSkeletonType& mesh_skeleton,
                            typename ProblemType::ProblemWriterType& writer) {
    using InterfaceTypeInternal = typename std::tuple_element<0, typename ProblemType::ProblemInterfaceTypes>::type;
    using InterfaceTypeLevee    = typename std::tuple_element<1, typename ProblemType::ProblemInterfaceTypes>::type;

    using EdgeInterfaceTypeInternal =
        typename std::tuple_element<0, typename ProblemType::ProblemEdgeInterfaceTypes>::type;
    using EdgeInterfaceTypeLevee =
        typename std::tuple_element<1, typename ProblemType::ProblemEdgeInterfaceTypes>::type;

    uint old_edge_internal = mesh_skeleton.GetNumberEdgeInterfaces();

    mesh.template CallForEachInterfaceOfType<InterfaceTypeInternal>([&mesh_skeleton](auto& intface) {
        mesh_skeleton.template CreateEdgeInterface<EdgeInterfaceTypeInternal>(intface);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of internal edges: "
                            << mesh_skeleton.GetNumberEdgeInterfaces() - old_edge_internal << std::endl;
    }

    uint old_edge_levee = mesh_skeleton.GetNumberEdgeInterfaces();

    mesh.template CallForEachInterfaceOfType<InterfaceTypeLevee>([&mesh_skeleton](auto& intface) {
        mesh_skeleton.template CreateEdgeInterface<EdgeInterfaceTypeLevee>(intface);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of levee edges: " << mesh_skeleton.GetNumberEdgeInterfaces() - old_edge_levee
                            << std::endl;
    }
}
}

#endif
