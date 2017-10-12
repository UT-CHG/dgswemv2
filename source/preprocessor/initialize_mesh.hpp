#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "geometry/mesh_definitions.hpp"

using namespace Geometry;

template <typename ProblemType>
void initialize_mesh_elements(typename ProblemType::ProblemMeshType& mesh, const MeshMetaData& mesh_data);

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh,
                                           Communicator& communicator,
                                           std::ofstream& log_file);

template <typename ProblemType, typename Communicator>
void initialize_mesh(typename ProblemType::ProblemMeshType& mesh,
                     const MeshMetaData& mesh_data,
                     Communicator& communicator,
                     const typename ProblemType::InputType& problem_specific_input,
                     std::ofstream& log_file) {
    initialize_mesh_elements<ProblemType>(mesh, mesh_data);

    log_file << "Number of elements: " << mesh.GetNumberElements() << std::endl;
    log_file << "Number of interfaces: " << mesh.GetNumberInterfaces() << std::endl;

    initialize_mesh_interfaces_boundaries<ProblemType, Communicator>(mesh, communicator, log_file);

    ProblemType::initialize_data_kernel(mesh, mesh_data, problem_specific_input);
}

template <typename ProblemType>
void initialize_mesh_elements(typename ProblemType::ProblemMeshType& mesh, const MeshMetaData& mesh_data) {
    using ElementType =
        typename std::tuple_element<0, Geometry::ElementTypeTuple<typename ProblemType::ProblemDataType>>::type;

    std::vector<Point<2>> nodal_coords_temp;
    for (const auto& element_meta : mesh_data.elements) {
        uint elt_id = element_meta.first;

        auto nodal_coordinates = mesh_data.GetNodalCoordinates(elt_id);

        for (auto& node_coordinate : nodal_coordinates) {
            nodal_coords_temp.push_back({node_coordinate[GlobalCoord::x], node_coordinate[GlobalCoord::y]});
        }

        mesh.template CreateElement<ElementType>(
            elt_id, nodal_coords_temp, element_meta.second.neighbor_ID, element_meta.second.boundary_type);
        nodal_coords_temp.clear();
    }
}

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh,
                                           Communicator& communicator,
                                           std::ofstream& log_file) {
    using RawBoundaryType = RawBoundary<1, typename ProblemType::ProblemDataType>;

    using InterfaceType =
        typename std::tuple_element<0, Geometry::InterfaceTypeTuple<typename ProblemType::ProblemDataType>>::type;

    std::map<uint, std::map<uint, RawBoundaryType>> pre_interfaces;
    std::map<uchar, std::vector<RawBoundaryType>> pre_boundaries;
    std::map<uint, std::map<uint, RawBoundaryType>> pre_distributed_boundaries;

    mesh.CallForEachElement([&pre_interfaces, &pre_boundaries, &pre_distributed_boundaries](auto& elem) {
        elem.CreateRawBoundaries(pre_interfaces, pre_boundaries, pre_distributed_boundaries);
    });

    for (auto it = pre_interfaces.begin(); it != pre_interfaces.end(); it++) {
        for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {

            mesh.template CreateInterface<InterfaceType>(itt->second, pre_interfaces.at(itt->first).at(it->first));

            pre_interfaces.at(itt->first).erase(it->first);
            it->second.erase(itt);
        }
    }

    ProblemType::create_boundaries_kernel(mesh, pre_boundaries, log_file);
    ProblemType::create_distributed_boundaries_kernel(mesh, communicator, pre_distributed_boundaries, log_file);
}
#endif
