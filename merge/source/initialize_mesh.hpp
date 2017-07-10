#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "geometry/mesh_definitions.hpp"

using namespace Geometry;

template<typename ProblemType>
void initialize_mesh_elements(typename ProblemType::mesh_type&);

template<typename ProblemType>
void initialize_mesh_interfaces_boundaries(typename ProblemType::mesh_type&);

template<typename ProblemType>
void initialize_mesh_VTK_geometry(typename ProblemType::mesh_type&);

template<typename ProblemType>
typename ProblemType::mesh_type* initialize_mesh(uint p) {
	typename ProblemType::mesh_type* mesh = new typename ProblemType::mesh_type(2);

	initialize_mesh_elements<ProblemType>(*mesh);
	initialize_mesh_interfaces_boundaries<ProblemType>(*mesh);
	initialize_mesh_VTK_geometry<ProblemType>(*mesh);

	return mesh;
}

template<typename ProblemType>
void initialize_mesh_elements(typename ProblemType::mesh_type& mesh) {
	using MasterType = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
	using ElementType = Element<2, MasterType, Shape::StraightTriangle, typename ProblemType::data_type>;

	double L = 90000;
	double W = 45000;

	uint m = 20;
	uint n = 10;

	double dx = L / m;
	double dy = W / n;

	uint ID;
	std::vector<Point<2>> nodal_coordinates(3);

	std::vector<unsigned char> boundaries(3);
	std::vector<uint> neighbors(3);

	for (uint i = 0; i < n; i++) {
		for (uint j = i % 2; j < m; j += 2) {
			ID = 2 * j + 2 * m * i;

			neighbors[0] = ID + 1;
			neighbors[1] = ID + 2 * m;
			neighbors[2] = ID - 1;

			boundaries[0] = SWE::internal;
			boundaries[1] = SWE::internal;
			boundaries[2] = SWE::internal;

			if (i == n - 1) {
				neighbors[1] = DEFAULT_ID;
				boundaries[1] = SWE::land;
			}
			if (j == 0) {
				neighbors[2] = DEFAULT_ID;
				boundaries[2] = SWE::land;
			}

			nodal_coordinates[0][0] = j*dx;
			nodal_coordinates[1][0] = nodal_coordinates[0][0];
			nodal_coordinates[2][0] = nodal_coordinates[0][0] + dx;

			nodal_coordinates[0][1] = (i + 1)*dy;
			nodal_coordinates[1][1] = nodal_coordinates[0][1] - dy;
			nodal_coordinates[2][1] = nodal_coordinates[0][1];

			mesh.template CreateElement<ElementType>(ID, nodal_coordinates, neighbors, boundaries);

			ID = ID + 1;

			neighbors[0] = ID - 1;
			neighbors[1] = ID + 1;
			neighbors[2] = ID - 2 * m;

			boundaries[0] = SWE::internal;
			boundaries[1] = SWE::internal;
			boundaries[2] = SWE::internal;

			if (i == 0) {
				neighbors[2] = DEFAULT_ID;
				boundaries[2] = SWE::land;
			}
			if (j == m - 1) {
				neighbors[1] = DEFAULT_ID;
				boundaries[1] = SWE::tidal;
			}

			nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

			nodal_coordinates[0][1] = nodal_coordinates[0][1] - dy;

			mesh.template CreateElement<ElementType>(ID, nodal_coordinates, neighbors, boundaries);
		}
	}

	for (uint i = 0; i < n; i++) {
		for (uint j = (i + 1) % 2; j < m; j += 2) {
			ID = 2 * j + 2 * m * i;

			neighbors[0] = ID + 1;
			neighbors[1] = ID - 1;
			neighbors[2] = ID - 2 * m;

			boundaries[0] = SWE::internal;
			boundaries[1] = SWE::internal;
			boundaries[2] = SWE::internal;

			if (i == 0) {
				neighbors[2] = DEFAULT_ID;
				boundaries[2] = SWE::land;
			}
			if (j == 0) {
				neighbors[1] = DEFAULT_ID;
				boundaries[1] = SWE::land;
			}

			nodal_coordinates[0][0] = j*dx;
			nodal_coordinates[1][0] = nodal_coordinates[0][0] + dx;
			nodal_coordinates[2][0] = nodal_coordinates[0][0];

			nodal_coordinates[0][1] = i*dy;
			nodal_coordinates[1][1] = nodal_coordinates[0][1];
			nodal_coordinates[2][1] = nodal_coordinates[0][1] + dy;

			mesh.template CreateElement<ElementType>(ID, nodal_coordinates, neighbors, boundaries);

			ID = ID + 1;

			neighbors[0] = ID - 1;
			neighbors[1] = ID + 2 * m;
			neighbors[2] = ID + 1;

			boundaries[0] = SWE::internal;
			boundaries[1] = SWE::internal;
			boundaries[2] = SWE::internal;

			if (i == n - 1) {
				neighbors[1] = DEFAULT_ID;
				boundaries[1] = SWE::land;
			}
			if (j == m - 1) {
				neighbors[2] = DEFAULT_ID;
				boundaries[2] = SWE::tidal;
			}

			nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

			nodal_coordinates[0][1] = nodal_coordinates[0][1] + dy;

			mesh.template CreateElement<ElementType>(ID, nodal_coordinates, neighbors, boundaries);
		}
	}
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

	ProblemType::create_boundaries(mesh, pre_boundaries);
}

template<typename ProblemType>
void initialize_mesh_VTK_geometry(typename ProblemType::mesh_type& mesh) {
	std::vector<Point<3>> points;
	Array2D<uint> cells;

	mesh.CallForEachElement(
		[&points, &cells](auto& elem)
	{ elem.InitializeVTK(points, cells); }
	);

	std::string file_name = "geometry.vtk";
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