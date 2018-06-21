#ifndef INITIALIZE_MESH_SKELETON_HPP
#define INITIALIZE_MESH_SKELETON_HPP

#include "input_parameters.hpp"
#include "geometry/mesh_definitions.hpp"
#include "simulation/writer.hpp"

template <typename ProblemType>
void initialize_mesh_skeleton(typename ProblemType::ProblemMeshType& mesh,
                              typename ProblemType::ProblemMeshSkeletonType& mesh_skeleton,
                              Writer<ProblemType>& writer) {
    using RawBoundaryType = Geometry::RawBoundary<1, typename ProblemType::ProblemDataType>;

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundaries;

    mesh.CallForEachElement([&raw_boundaries](auto& elem) { elem.CreateRawBoundaries(raw_boundaries); });

    using EdgeInternalType =
        typename std::tuple_element<0,
                                    Geometry::EdgeInternalTypeTuple<typename ProblemType::ProblemDataType,
                                                                    typename ProblemType::ProblemEdgeDataType>>::type;

    using EdgeBoundaryType =
        typename std::tuple_element<0,
                                    Geometry::EdgeBoundaryTypeTuple<typename ProblemType::ProblemDataType,
                                                                    typename ProblemType::ProblemEdgeDataType>>::type;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        if (is_internal(it->first)) {
            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                if (it->second.find(key_pre_int_ex) != it->second.end()) {
                    auto& raw_boundary_in = itt->second;
                    auto& raw_boundary_ex = it->second.find(key_pre_int_ex)->second;

                    mesh_skeleton.template CreateEdgeInternal<EdgeInternalType>(raw_boundary_in, raw_boundary_ex);
                }

                it->second.erase(itt++);
            }
        } else {
            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                auto& raw_boundary = itt->second;

                mesh_skeleton.template CreateEdgeBoundary<EdgeBoundaryType>(raw_boundary);

                it->second.erase(itt++);
            }
        }
    }

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        if (it->second.size() != 0) {
            throw std::logic_error("Fatal Error: unprocessed raw_boundaries of boundary type " +
                                   std::to_string(it->first) + "!\n");
        }
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of edges internal: " << mesh_skeleton.GetNumberEdgeInternals() << std::endl;
        writer.GetLogFile() << "Number of edges boundary: " << mesh_skeleton.GetNumberEdgeBoundaries() << std::endl;
    }
}

#endif