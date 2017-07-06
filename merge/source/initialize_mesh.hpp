#ifndef INITIALIZE_MESH_HPP
#define INITIALIZE_MESH_HPP

#include "general_definitions.hpp"

#include "mesh_definitions.hpp"

namespace Geometry {	
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

				//mesh.template CreateElement<ElementType>(ID, *triangle, nodal_coordinates, neighbors, boundaries);

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

				//this->elements[ID] = new Element<>(ID, *triangle, nodal_coordinates, neighbors, boundaries);
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

				//this->elements[ID] = new Element<>(ID, *triangle, nodal_coordinates, neighbors, boundaries);

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

				//this->elements[ID] = new Element<>(ID, *triangle, nodal_coordinates, neighbors, boundaries);
			}
		}



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
}

#endif