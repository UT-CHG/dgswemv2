#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "general_definitions.hpp"

#include "mesh_definitions.hpp"

namespace Geometry {
	template<typename Data>
	void initialize_mesh_VTK_geometry(MeshType<Data>&);
	
	template<typename Data>
	void initialize_mesh_interfaces_boundaries(MeshType<Data>&);

	template<typename Data>
	void initialize_mesh(int p, MeshType<Data>& mesh) {
		using MasterType = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
		using ElementType = Element<2, MasterType, Shape::StraightTriangle, Data>;

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

		MasterType* triangle = new MasterType(p);
		
		for (uint i = 0; i < n; i++) {
			for (uint j = i % 2; j < m; j += 2) {
				ID = 2 * j + 2 * m * i;

				neighbors[0] = ID + 1;
				neighbors[1] = ID + 2 * m;
				neighbors[2] = ID - 1;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == n - 1) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}
				if (j == 0) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}

				nodal_coordinates[0][0] = j*dx;
				nodal_coordinates[1][0] = nodal_coordinates[0][0];
				nodal_coordinates[2][0] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[0][1] = (i + 1)*dy;
				nodal_coordinates[1][1] = nodal_coordinates[0][1] - dy;
				nodal_coordinates[2][1] = nodal_coordinates[0][1];

				mesh.template CreateElement<ElementType>(ID, ID, *triangle, nodal_coordinates, neighbors, boundaries);

				ID = ID + 1;

				neighbors[0] = ID - 1;
				neighbors[1] = ID + 1;
				neighbors[2] = ID - 2 * m;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == 0) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}
				if (j == m - 1) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = OCEAN;
				}

				nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[0][1] = nodal_coordinates[0][1] - dy;

				mesh.template CreateElement<ElementType>(ID, ID, *triangle, nodal_coordinates, neighbors, boundaries);
			}
		}

		for (uint i = 0; i < n; i++) {
			for (uint j = (i + 1) % 2; j < m; j += 2) {
				ID = 2 * j + 2 * m * i;

				neighbors[0] = ID + 1;
				neighbors[1] = ID - 1;
				neighbors[2] = ID - 2 * m;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == 0) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = LAND;
				}
				if (j == 0) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}

				nodal_coordinates[0][0] = j*dx;
				nodal_coordinates[1][0] = nodal_coordinates[0][0] + dx;
				nodal_coordinates[2][0] = nodal_coordinates[0][0];

				nodal_coordinates[0][1] = i*dy;
				nodal_coordinates[1][1] = nodal_coordinates[0][1];
				nodal_coordinates[2][1] = nodal_coordinates[0][1] + dy;

				mesh.template CreateElement<ElementType>(ID, ID, *triangle, nodal_coordinates, neighbors, boundaries);

				ID = ID + 1;

				neighbors[0] = ID - 1;
				neighbors[1] = ID + 2 * m;
				neighbors[2] = ID + 1;

				boundaries[0] = INTERNAL;
				boundaries[1] = INTERNAL;
				boundaries[2] = INTERNAL;

				if (i == n - 1) {
					neighbors[1] = DEFAULT_ID;
					boundaries[1] = LAND;
				}
				if (j == m - 1) {
					neighbors[2] = DEFAULT_ID;
					boundaries[2] = OCEAN;
				}

				nodal_coordinates[0][0] = nodal_coordinates[0][0] + dx;

				nodal_coordinates[0][1] = nodal_coordinates[0][1] + dy;

				mesh.template CreateElement<ElementType>(ID, ID, *triangle, nodal_coordinates, neighbors, boundaries);
			}
		}


		initialize_mesh_VTK_geometry<Data>(mesh);
		initialize_mesh_interfaces_boundaries<Data>(mesh);

		printf("%d\n", mesh.GetNumberElements());
		printf("%d\n", mesh.GetNumberInterfaces());

		//mesh.CallForEachElement([](auto& elem){ printf("%d\n", elem.ID);});
		

		/*
		using ElementType = Element::Triangle< Basis::Dubiner, typename Data::Volume>;
		{//make all the elements
			using Point = std::array<double, 2>;

			using MasterTriangle = Element::MasterTriangle<Basis::Dubiner>;
			std::shared_ptr<MasterTriangle > master(std::make_shared<MasterTriangle>(p));


			for (const auto& elt : mesh_file.elements) {
				std::vector<Point> vrtxs(3);
				for (uint i = 0; i < 3; ++i) {
					std::array<double, 3> tmp = mesh_file.nodes.at(elt.second[i + 1]);
					vrtxs[i][0] = tmp[0];
					vrtxs[i][1] = tmp[1];
				}

				Shape::Straight<2> shape(vrtxs);
				mesh.template create_element<ElementType>(elt.first, master, shape);
			}
		}

		{//make all edges
			using eltID_faceID = std::pair<int, int>;

			std::unordered_map<std::uint64_t, std::pair<eltID_faceID, eltID_faceID> > edge_dictionary;
			for (const auto& elt : mesh_file.elements) {
				std::vector<int> node{ elt.second[1], elt.second[2], elt.second[3] };

				for (int k = 0; k < 3; ++k) {
					std::uint64_t curr_key = (static_cast<std::uint64_t>(std::min(node[k], node[(k + 1) % 3]))) << 32
						| static_cast<std::uint64_t>(std::max(node[k], node[(k + 1) % 3]));

					if (edge_dictionary.count(curr_key)) {//if already one element has the edge.
						edge_dictionary.at(curr_key).second = std::make_pair(elt.first, k);
					}
					else {
						std::pair<eltID_faceID, eltID_faceID> edge_info{ {elt.first,k},{-1,0} };
						edge_dictionary.insert({ curr_key, edge_info });
					}
				}
			}

			for (const auto& edge : edge_dictionary) {
				//check if there are two elements associated with this edge
				if (edge.second.second.first != -1) {
					int eltA_id = edge.second.first.first;
					int eltB_id = edge.second.second.first;
					int fidA = edge.second.first.second;
					int fidB = edge.second.second.second;

					mesh.template create_interior_edge<Edge::InteriorEdge,
						typename Data::Edge,
						ElementType, ElementType>(eltA_id, eltB_id,
							fidA, fidB);
				}
				else {
					//treat boundary conditions
					int elt_id = edge.second.first.first;
					int fid = edge.second.first.second;

					mesh.template create_boundary_edge<Edge::BoundaryEdge,
						BoundaryData,
						ElementType>(elt_id, fid);
				}
			}
		}
		*/
	}


template<typename Data>
void initialize_mesh_interfaces_boundaries(MeshType<Data>& mesh) {
	using RawBoundaryType = RawBoundary<1>;

	using InterfaceType = Interface<1, Integration::GaussLegendre_1D, Data>;

	std::map<uint, std::map<uint, RawBoundaryType*>> pre_interfaces;
	std::map<unsigned char, std::vector<RawBoundaryType*>> pre_boundaries;
	
	mesh.CallForEachElement([&pre_interfaces, &pre_boundaries](auto& elem){ elem.CreateRawBoundaries(pre_interfaces, pre_boundaries);});

/*
	for (auto it = pre_boundaries.begin(); it != pre_boundaries.end(); it++) {
		this->boundaries[it->first] = std::vector<Boundary<>*>();

		for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
			this->boundaries[it->first].push_back(new Boundary<>(*(*itt)));
			delete *itt;
		}
	}

	*/
	
	for (auto it = pre_interfaces.begin(); it != pre_interfaces.end(); it++) {
		for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {

			mesh.template CreateInterface<InterfaceType>(
				*itt->second,
				*pre_interfaces[itt->first][it->first]);

			delete itt->second;
			delete pre_interfaces[itt->first][it->first];

			pre_interfaces[itt->first].erase(it->first);
			it->second.erase(itt);
		}
	}
}

template<typename Data>
void initialize_mesh_VTK_geometry(MeshType<Data>& mesh) {
    std::vector<Point<3>> points;
    Array2D<uint> cells;

	mesh.CallForEachElement([&points, &cells](auto& elem){ elem.InitializeVTK(points,cells);});

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
}

#endif