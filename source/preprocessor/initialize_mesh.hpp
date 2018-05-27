#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "input_parameters.hpp"
#include "../geometry/mesh_definitions.hpp"
#include "../simulation/writer.hpp"

template <typename ProblemType>
void initialize_mesh_elements(typename ProblemType::ProblemMeshType& mesh,
                              InputParameters<typename ProblemType::ProblemInputType>& input,
                              Writer<ProblemType>& writer);

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh,
                                           InputParameters<typename ProblemType::ProblemInputType>& input,
                                           Communicator& communicator,
                                           Writer<ProblemType>& writer);

template <typename ProblemType, typename Communicator>
void initialize_mesh(typename ProblemType::ProblemMeshType& mesh,
                     InputParameters<typename ProblemType::ProblemInputType>& input,
                     Communicator& communicator,
                     Writer<ProblemType>& writer) {
    mesh.SetMeshName(input.mesh_input.mesh_data.mesh_name);

    initialize_mesh_elements<ProblemType>(mesh, input, writer);

    initialize_mesh_interfaces_boundaries<ProblemType, Communicator>(mesh, input, communicator, writer);
}

template <typename ProblemType>
void initialize_mesh_elements(typename ProblemType::ProblemMeshType& mesh,
                              InputParameters<typename ProblemType::ProblemInputType>& input,
                              Writer<ProblemType>& writer) {
    MeshMetaData& mesh_data = input.mesh_input.mesh_data;

    using ElementType =
        typename std::tuple_element<0, Geometry::ElementTypeTuple<typename ProblemType::ProblemDataType>>::type;

    std::vector<Point<2>> nodal_coords_temp;
    for (auto& element_meta : mesh_data.elements) {
        uint elt_id = element_meta.first;

        auto nodal_coordinates = mesh_data.get_nodal_coordinates(elt_id);

        for (auto& node_coordinate : nodal_coordinates) {
            nodal_coords_temp.push_back({node_coordinate[GlobalCoord::x], node_coordinate[GlobalCoord::y]});
        }

        mesh.template CreateElement<ElementType>(elt_id,
                                                 nodal_coords_temp,
                                                 element_meta.second.node_ID,
                                                 element_meta.second.neighbor_ID,
                                                 element_meta.second.boundary_type);
        nodal_coords_temp.clear();
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of elements: " << mesh.GetNumberElements() << std::endl;
    }
}

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh,
                                           InputParameters<typename ProblemType::ProblemInputType>& input,
                                           Communicator& communicator,
                                           Writer<ProblemType>& writer) {
    using RawBoundaryType = Geometry::RawBoundary<1, typename ProblemType::ProblemDataType>;

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundaries;

    mesh.CallForEachElement([&raw_boundaries](auto& elem) { elem.CreateRawBoundaries(raw_boundaries); });

    ProblemType::create_interfaces_kernel(raw_boundaries, mesh, input, writer);
    ProblemType::create_boundaries_kernel(raw_boundaries, mesh, input, writer);
    ProblemType::create_distributed_boundaries_kernel(raw_boundaries, mesh, input, communicator, writer);

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        if (it->second.size() != 0) {
            throw std::logic_error("Error: Unprocessed raw_boundaries of boundary type " + std::to_string(it->first) +
                                   "\n");
        }
    }
}

#endif