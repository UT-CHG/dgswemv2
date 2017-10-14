#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "geometry/mesh_definitions.hpp"

using namespace Geometry;

template <typename ProblemType>
void initialize_mesh_elements(typename ProblemType::ProblemMeshType& mesh, const MeshMetaData& mesh_data);

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh, Communicator& communicator);

template <typename ProblemType>
void initialize_mesh_VTK_geometry(typename ProblemType::ProblemMeshType& mesh);

template <typename ProblemType, typename Communicator>
void initialize_mesh(typename ProblemType::ProblemMeshType& mesh,
                     const MeshMetaData& mesh_data,
                     Communicator& communicator) {
    initialize_mesh_elements<ProblemType>(mesh, mesh_data);
    initialize_mesh_interfaces_boundaries<ProblemType, Communicator>(mesh, communicator);
#ifdef OUTPUT
    initialize_mesh_VTK_geometry<ProblemType>(mesh);
#endif
    ProblemType::initialize_data_kernel(mesh, mesh_data);
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
#ifdef VERBOSE
    std::ofstream log_file("output/" + mesh.GetMeshName() + "_log", std::ofstream::app);

    log_file << "Number of elements: " << mesh.GetNumberElements() << std::endl;
#endif
}

template <typename ProblemType, typename Communicator>
void initialize_mesh_interfaces_boundaries(typename ProblemType::ProblemMeshType& mesh, Communicator& communicator) {
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
#ifdef VERBOSE
    std::ofstream log_file("output/" + mesh.GetMeshName() + "_log", std::ofstream::app);

    log_file << "Number of interfaces: " << mesh.GetNumberInterfaces() << std::endl;
#endif
    ProblemType::create_boundaries_kernel(mesh, pre_boundaries);
    ProblemType::create_distributed_boundaries_kernel(mesh, communicator, pre_distributed_boundaries);
}

template <typename ProblemType>
void initialize_mesh_VTK_geometry(typename ProblemType::ProblemMeshType& mesh) {
    std::vector<Point<3>> points;
    Array2D<uint> cells;

    mesh.CallForEachElement([&points, &cells](auto& elem) { elem.InitializeVTK(points, cells); });

    std::string file_name = "output/" + mesh.GetMeshName() + "_geometry.vtk";
    std::ofstream file(file_name);

    file << "# vtk DataFile Version 3.0\n";
    file << "OUTPUT DATA\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << points.size() << " float\n";

    for (auto it = points.begin(); it != points.end(); it++) {
        file << (*it)[0] << '\t' << (*it)[1] << '\t' << (*it)[2] << '\n';
    }

    uint n_cell_entries = 0;
    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
            case VTKElementTypes::straight_triangle:
                n_cell_entries += 4;
                break;
            default:
                printf("\n");
                printf("MESH InitializeVTK - Fatal error!\n");
                printf("Undefined cell type = %d\n", (*it)[0]);
                exit(1);
        }
    }

    file << "CELLS " << cells.size() << ' ' << n_cell_entries << '\n';

    uint n_nodes;

    for (auto it = cells.begin(); it != cells.end(); it++) {
        switch ((*it)[0]) {
            case VTKElementTypes::straight_triangle:
                file << 3 << '\t';
                n_nodes = 3;
                break;
            default:
                printf("\n");
                printf("MESH InitializeVTK - Fatal error!\n");
                printf("Undefined cell type = %d\n", (*it)[0]);
                exit(1);
        }

        for (uint i = 1; i <= n_nodes; i++) {
            file << (*it)[i] << '\t';
        }
        file << '\n';
    }

    file << "CELL_TYPES " << cells.size() << '\n';

    for (auto it = cells.begin(); it != cells.end(); it++) {
        file << (*it)[0] << '\n';
    }
}

#endif
