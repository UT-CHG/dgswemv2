#ifndef MESH_METADATA_HPP
#define MESH_METADATA_HPP

#include "general_definitions.hpp"

#include "ADCIRC_reader/adcirc_format.hpp"
#include "problem/SWE/swe_definitions.hpp"

struct NodeMetaData {
        Point<2> coordinates;
        double bathymetry;
};


// ID: corresponds to the element ID
// coordinates: correspond to the vertices of the element
//   moving counter clockwise around the triangle
// boundary type: coresponds to the type of the element edge (Interface, Tidal boundary, etc.)
// neighbor_IDs: corresponds to the neighbor ID if it exists. If it does not exist the default ID is set to
//   DEFAULT_ID as set in general_definitions.hpp
struct ElementMetaData {
	ElementMetaData(uint n_faces)
		: node_ids(n_faces), neighbor_ID(n_faces), boundary_type(n_faces)
	{}

	std::vector<uint> node_ids;
	std::vector<uint> neighbor_ID;
	std::vector<uchar> boundary_type;
};

struct MeshMetaData {
        MeshMetaData() = default;
	MeshMetaData(const AdcircFormat& mesh_file);\

        std::vector<Point<2> > GetNodalCoordinates(uint elt_id) const;
        std::vector<double> GetBathymetry(uint elt_id) const;

	std::unordered_map<uint, ElementMetaData> _elements;
        std::unordered_map<uint, NodeMetaData> _nodes;
};

#endif
