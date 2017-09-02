#ifndef MESH_METADATA_HPP
#define MESH_METADATA_HPP

#include "../general_definitions.hpp"

#include "ADCIRC_reader/adcirc_format.hpp"
#include "../problem/SWE/swe_definitions.hpp"

struct NodeMetaData {
    Point<3> coordinates;

    friend std::ostream& operator<<(std::ostream& s, const NodeMetaData& node) {
        return s << std::setprecision(15) << node.coordinates[0] << " " << node.coordinates[1] << " "
                 << node.coordinates[3];
    }

    friend std::istream& operator>>(std::istream& s, NodeMetaData& node) {
        return s >> node.coordinates[0] >> node.coordinates[1] >> node.coordinates[3];
    }

    friend bool operator==(const NodeMetaData& lhs, const NodeMetaData& rhs) {
        return (lhs.coordinates[0] == rhs.coordinates[0]) && (lhs.coordinates[1] == rhs.coordinates[1]) &&
               (lhs.coordinates[3] == rhs.coordinates[3]);
    }
};

// ID: corresponds to the element ID
// coordinates: correspond to the vertices of the element
//   moving counter clockwise around the triangle
// boundary type: coresponds to the type of the element edge (Interface, Tidal boundary, etc.)
// neighbor_IDs: corresponds to the neighbor ID if it exists. If it does not exist the default ID is set to
//   DEFAULT_ID as set in general_definitions.hpp
struct ElementMetaData {
    ElementMetaData() = default;
    ElementMetaData(uint n_faces) : node_ids(n_faces), neighbor_ID(n_faces), boundary_type(n_faces) {}

    std::vector<uint> node_ids;
    std::vector<uint> neighbor_ID;
    std::vector<uchar> boundary_type;
    friend std::ostream& operator<<(std::ostream& s, const ElementMetaData& elt) {
        s << elt.node_ids.size();
        for (const auto& node_id : elt.node_ids) {
            s << " " << node_id;
        }

        for (const auto& neigh_id : elt.neighbor_ID) {
            s << " " << neigh_id;
        }

        for (const auto& bt : elt.boundary_type) {
            s << " " << bt;
        }

        return s;
    }

    friend std::istream& operator>>(std::istream& s, ElementMetaData& elt) {
        uint n_faces;
        s >> n_faces;

        elt.node_ids.resize(n_faces);
        elt.neighbor_ID.resize(n_faces);
        elt.boundary_type.resize(n_faces);

        for (uint i = 0; i < n_faces; ++i) {
            s >> elt.node_ids[i];
        }

        for (uint i = 0; i < n_faces; ++i) {
            s >> elt.neighbor_ID[i];
        }

        for (uint i = 0; i < n_faces; ++i) {
            s >> elt.boundary_type[i];
        }

        return s;
    }
};

inline bool operator==(const ElementMetaData& lhs, const ElementMetaData& rhs) {
    return (lhs.node_ids == rhs.node_ids) && (lhs.neighbor_ID == rhs.neighbor_ID) &&
           (lhs.boundary_type == rhs.boundary_type);
}

struct MeshMetaData {
    MeshMetaData() = default;  // why define default constructor?
    MeshMetaData(const AdcircFormat& mesh_file);
    MeshMetaData(const std::string& file);  // read from file

    void WriteTo(const std::string& file);  // write to file

    std::vector<Point<3>> GetNodalCoordinates(uint elt_id) const;

    std::string _mesh_name;
    std::unordered_map<uint, ElementMetaData> _elements;
    std::unordered_map<uint, NodeMetaData> _nodes;
};

struct DistributedBoundaryMetaData {
    std::pair<uint, uint> elements;
    std::pair<uint, uint> bound_ids;
    uint p;

    friend std::ostream& operator<<(std::ostream& s, const DistributedBoundaryMetaData& dist_int) {
        return s << dist_int.elements.first << " " << dist_int.elements.second << " " << dist_int.bound_ids.first << " "
                 << dist_int.bound_ids.second << " " << dist_int.p;
    }

    friend std::istream& operator>>(std::istream& s, DistributedBoundaryMetaData& dist_int) {
        return s >> dist_int.elements.first >> dist_int.elements.second >> dist_int.bound_ids.first >>
               dist_int.bound_ids.second >> dist_int.p;
    }
};

#endif
