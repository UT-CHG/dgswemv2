#ifndef MESH_METADATA_HPP
#define MESH_METADATA_HPP

#include "general_definitions.hpp"
#include "ADCIRC_reader/adcirc_format.hpp"
#include "utilities/file_exists.hpp"

struct NodeMetaData {
    Point<3> coordinates;

    friend std::ostream& operator<<(std::ostream& s, const NodeMetaData& node) {
        return s << std::setprecision(15) << node.coordinates[0] << " " << node.coordinates[1] << " "
                 << node.coordinates[2];
    }

    friend std::istream& operator>>(std::istream& s, NodeMetaData& node) {
        return s >> node.coordinates[0] >> node.coordinates[1] >> node.coordinates[2];
    }

    friend bool operator==(const NodeMetaData& lhs, const NodeMetaData& rhs) {
        return (lhs.coordinates[0] == rhs.coordinates[0]) && (lhs.coordinates[1] == rhs.coordinates[1]) &&
               (lhs.coordinates[2] == rhs.coordinates[2]);
    }
};

struct ElementMetaData {
    ElementMetaData() = default;
    ElementMetaData(uint n_faces) : node_ID(n_faces), neighbor_ID(n_faces), boundary_type(n_faces) {}

    std::vector<uint> node_ID;
    std::vector<uint> neighbor_ID;
    std::vector<uchar> boundary_type;

    friend std::ostream& operator<<(std::ostream& s, const ElementMetaData& elt) {
        s << elt.node_ID.size();
        for (const auto& node_id : elt.node_ID) {
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

        elt.node_ID.resize(n_faces);
        elt.neighbor_ID.resize(n_faces);
        elt.boundary_type.resize(n_faces);

        for (uint i = 0; i < n_faces; ++i) {
            s >> elt.node_ID[i];
        }

        for (uint i = 0; i < n_faces; ++i) {
            s >> elt.neighbor_ID[i];
        }

        for (uint i = 0; i < n_faces; ++i) {
            s >> elt.boundary_type[i];
        }

        return s;
    }

    friend bool operator==(const ElementMetaData& lhs, const ElementMetaData& rhs) {
        return (lhs.node_ID == rhs.node_ID) && (lhs.neighbor_ID == rhs.neighbor_ID) &&
               (lhs.boundary_type == rhs.boundary_type);
    }
};

struct MeshMetaData {
    MeshMetaData() = default;
    MeshMetaData(const AdcircFormat& mesh_file);
    MeshMetaData(const std::string& mesh_file);  // read from file

    void write_to(const std::string& file);  // write to file

    std::vector<Point<3>> get_nodal_coordinates(uint elt_id) const;

    std::string mesh_name;
    std::unordered_map<uint, ElementMetaData> elements;
    std::unordered_map<uint, NodeMetaData> nodes;
};

struct DBPairMetaData {
    std::pair<uint, uint> elements;
    std::pair<uint, uint> bound_ids;
    uint p;

    friend std::ostream& operator<<(std::ostream& s, const DBPairMetaData& dist_int) {
        return s << dist_int.elements.first << " " << dist_int.elements.second << " " << dist_int.bound_ids.first << " "
                 << dist_int.bound_ids.second << " " << dist_int.p;
    }

    friend std::istream& operator>>(std::istream& s, DBPairMetaData& dist_int) {
        return s >> dist_int.elements.first >> dist_int.elements.second >> dist_int.bound_ids.first >>
               dist_int.bound_ids.second >> dist_int.p;
    }
};

struct RankBoundaryMetaData {
    uint locality_in;
    uint locality_ex;

    uint submesh_in;
    uint submesh_ex;

    std::vector<uint> elements_in;
    std::vector<uint> elements_ex;

    std::vector<uint> bound_ids_in;
    std::vector<uint> bound_ids_ex;

    std::vector<uint> p;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & locality_in
            & locality_ex
            & submesh_in
            & submesh_ex
            & elements_in
            & elements_ex
            & bound_ids_in
            & bound_ids_ex
            & p;
        // clang-format on
    }
#endif
};

struct DistributedBoundaryMetaData {
    std::vector<RankBoundaryMetaData> rank_boundary_data;

    DistributedBoundaryMetaData() = default;
    DistributedBoundaryMetaData(const std::string& dbmd_file, uint locality_id, uint submesh_id);  // read from file
};

#endif