#ifndef INITIALIZE_MESH_SKELETON_HPP
#define INITIALIZE_MESH_SKELETON_HPP

#include "input_parameters.hpp"
#include "geometry/mesh_definitions.hpp"
#include "simulation/writer.hpp"

template <typename ProblemType>
void initialize_mesh_skeleton(typename ProblemType::ProblemMeshType& mesh,
                              typename ProblemType::ProblemMeshSkeletonType& mesh_skeleton,
                              Writer<ProblemType>& writer) {
    using EdgeInternalType =
        typename std::tuple_element<0,
                                    Geometry::EdgeInternalTypeTuple<typename ProblemType::ProblemDataType,
                                                                    typename ProblemType::ProblemEdgeDataType>>::type;

    using EdgeBoundaryType =
        typename std::tuple_element<0,
                                    Geometry::EdgeBoundaryTypeTuple<typename ProblemType::ProblemDataType,
                                                                    typename ProblemType::ProblemEdgeDataType>>::type;

    mesh.CallForEachInterface(
        [&mesh_skeleton](auto& intface) { mesh_skeleton.template CreateEdgeInternal<EdgeInternalType>(intface); });

    mesh.CallForEachBoundary(
        [&mesh_skeleton](auto& bound) { mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryType>(bound); });

    mesh.CallForEachDistributedBoundary(
        [&mesh_skeleton](auto& dbound) { mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryType>(dbound); });

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of edges internal: " << mesh_skeleton.GetNumberEdgeInternals() << std::endl;
        writer.GetLogFile() << "Number of edges boundary: " << mesh_skeleton.GetNumberEdgeBoundaries() << std::endl;
    }
}

#endif