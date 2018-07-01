#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "input_parameters.hpp"
#include "geometry/mesh_definitions.hpp"
#include "simulation/writer.hpp"

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

    initialize_mesh_interfaces_boundaries<ProblemType, Communicator>(mesh, input.problem_input, communicator, writer);
}

template <typename ProblemType>
void initialize_mesh_elements(typename ProblemType::ProblemMeshType& mesh,
                              InputParameters<typename ProblemType::ProblemInputType>& input,
                              Writer<ProblemType>& writer) {
    MeshMetaData& mesh_data = input.mesh_input.mesh_data;

    using ElementType =
        typename std::tuple_element<0, Geometry::ElementTypeTuple<typename ProblemType::ProblemDataType>>::type;

    for (auto& element_meta : mesh_data.elements) {
        uint elt_id = element_meta.first;

        auto nodal_coordinates = mesh_data.get_nodal_coordinates(elt_id);

        mesh.template CreateElement<ElementType>(elt_id,
                                                 std::move(nodal_coordinates),
                                                 std::move(element_meta.second.node_ID),
                                                 std::move(element_meta.second.neighbor_ID),
                                                 std::move(element_meta.second.boundary_type));
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of elements: " << mesh.GetNumberElements() << std::endl;
    }
}

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh,
                                           typename ProblemType::ProblemInputType& problem_input,
                                           Communicator& communicator,
                                           Writer<ProblemType>& writer) {
    using RawBoundaryType = Geometry::RawBoundary<1, typename ProblemType::ProblemDataType>;

    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>> raw_boundaries;

    mesh.CallForEachElement([&raw_boundaries](auto& elem) { elem.CreateRawBoundaries(raw_boundaries); });

    ProblemType::create_interfaces_kernel(raw_boundaries, mesh, problem_input, writer);
    ProblemType::create_boundaries_kernel(raw_boundaries, mesh, problem_input, writer);
    ProblemType::create_distributed_boundaries_kernel(raw_boundaries, mesh, problem_input, communicator, writer);

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        if (it->second.size() != 0) {
            throw std::logic_error("Fatal Error: unprocessed raw_boundaries of boundary type " +
                                   std::to_string(it->first) + "!\n");
        }
    }
}

#endif