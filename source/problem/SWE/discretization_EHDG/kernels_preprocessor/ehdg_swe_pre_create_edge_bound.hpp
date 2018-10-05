#ifndef EHDG_SWE_PRE_CREATE_EDGE_BOUND_HPP
#define EHDG_SWE_PRE_CREATE_EDGE_BOUND_HPP

namespace SWE {
namespace EHDG {
void Problem::create_edge_boundaries(ProblemMeshType& mesh,
                                     ProblemMeshSkeletonType& mesh_skeleton,
                                     ProblemWriterType& writer) {
    using BoundaryTypes = Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow>;

    using EdgeBoundaryTypes =
        Geometry::EdgeBoundaryTypeTuple<EdgeData,
                                        Geometry::BoundaryTypeTuple<Data, BC::Land, BC::Tide, BC::Flow>>::Type;

    uint old_edge_land = mesh_skeleton.GetNumberEdgeBoundaries();

    using BoundaryTypeLand = typename std::tuple_element<0, BoundaryTypes>::type;
    mesh.CallForEachBoundaryOfType<BoundaryTypeLand>([&mesh_skeleton, &writer](auto& bound) {
        using EdgeBoundaryTypeLand = typename std::tuple_element<0, EdgeBoundaryTypes>::type;

        mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryTypeLand>(bound);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of land edges: " << mesh_skeleton.GetNumberEdgeBoundaries() - old_edge_land
                            << std::endl;
    }

    uint old_edge_tide = mesh_skeleton.GetNumberEdgeBoundaries();

    using BoundaryTypeTide = typename std::tuple_element<1, BoundaryTypes>::type;
    mesh.CallForEachBoundaryOfType<BoundaryTypeTide>([&mesh_skeleton, &writer](auto& bound) {
        using EdgeBoundaryTypeTide = typename std::tuple_element<1, EdgeBoundaryTypes>::type;

        mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryTypeTide>(bound);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of tide edges: " << mesh_skeleton.GetNumberEdgeBoundaries() - old_edge_tide
                            << std::endl;
    }

    uint old_edge_flow = mesh_skeleton.GetNumberEdgeBoundaries();

    using BoundaryTypeFlow = typename std::tuple_element<2, BoundaryTypes>::type;
    mesh.CallForEachBoundaryOfType<BoundaryTypeFlow>([&mesh_skeleton, &writer](auto& bound) {
        using EdgeBoundaryTypeFlow = typename std::tuple_element<2, EdgeBoundaryTypes>::type;

        mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryTypeFlow>(bound);
    });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of flow edges: " << mesh_skeleton.GetNumberEdgeBoundaries() - old_edge_flow
                            << std::endl;
    }
}
}
}

#endif
