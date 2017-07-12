#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/swe_definitions.hpp"

using namespace Geometry;

template<typename ProblemType>
void initialize_mesh_elements(typename ProblemType::mesh_type&, const MeshMetaData& mesh_data);

template<typename ProblemType>
void initialize_mesh_interfaces_boundaries(typename ProblemType::mesh_type&);

template<typename ProblemType>
void initialize_mesh_VTK_geometry(typename ProblemType::mesh_type&);

template<typename ProblemType>
typename ProblemType::mesh_type* initialize_mesh(uint p, const MeshMetaData& mesh_data) {
	typename ProblemType::mesh_type* mesh = new typename ProblemType::mesh_type(2);

	initialize_mesh_elements<ProblemType>(*mesh, mesh_data);
	initialize_mesh_interfaces_boundaries<ProblemType>(*mesh);
	initialize_mesh_VTK_geometry<ProblemType>(*mesh);

	return mesh;
}

template<typename ProblemType>
void initialize_mesh_elements(typename ProblemType::mesh_type& mesh, const MeshMetaData& mesh_data) {
	using MasterType = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
	using ElementType = Element<2, MasterType, Shape::StraightTriangle, typename ProblemType::data_type>;

	for (auto it = mesh_data._meta.begin(); it != mesh_data._meta.end(); it++) {
		mesh.template CreateElement<ElementType>(it->first, it->second.nodal_coordinates,
			it->second.neighbor_ID, it->second.boundary_type);
	}

	std::cout << "Number of elements: " << mesh.GetNumberElements() << "\n";
}

template<typename ProblemType>
void initialize_mesh_interfaces_boundaries(typename ProblemType::mesh_type& mesh) {
	using RawBoundaryType = RawBoundary<1, typename ProblemType::data_type>;

	using InterfaceType = Interface<1, Integration::GaussLegendre_1D, typename ProblemType::data_type>;

	std::map<uint, std::map<uint, RawBoundaryType>> pre_interfaces;
	std::map<unsigned char, std::vector<RawBoundaryType>> pre_boundaries;

	mesh.CallForEachElement(
		[&pre_interfaces, &pre_boundaries](auto& elem)
		{ elem.CreateRawBoundaries(pre_interfaces, pre_boundaries); }
	);

	for (auto it = pre_interfaces.begin(); it != pre_interfaces.end(); it++) {
		for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {

			mesh.template CreateInterface<InterfaceType>(
				itt->second,
				pre_interfaces.at(itt->first).at(it->first)
				);

			pre_interfaces.at(itt->first).erase(it->first);
			it->second.erase(itt);
		}
	}

	std::cout << "Number of interfaces: " << mesh.GetNumberInterfaces() << "\n";

	ProblemType::create_boundaries_kernel(mesh, pre_boundaries);
}

template<typename ProblemType>
void initialize_mesh_VTK_geometry(typename ProblemType::mesh_type& mesh) {
	std::vector<Point<3>> points;
	Array2D<uint> cells;

	mesh.CallForEachElement(
		[&points, &cells](auto& elem)
		{ elem.InitializeVTK(points, cells); }
	);

	std::string file_name = "output/geometry.vtk";
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
		case TRIANGLE: n_cell_entries += 4; break;
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
		case TRIANGLE: file << 3 << '\t'; n_nodes = 3; break;
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