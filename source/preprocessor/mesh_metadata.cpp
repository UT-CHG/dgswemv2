#include "mesh_metadata.hpp"

MeshMetaData::MeshMetaData(const AdcircFormat& mesh_file) {
    for (const auto& nod : mesh_file.nodes) {
        if (nod.first < 0) {
            throw std::logic_error("ERROR in mesh_metadata.cpp: Node ID is negative; not supported\n");
        }
        uint ID = nod.first;
        _nodes.insert({ID, {{nod.second[0], nod.second[1]}, nod.second[2]}});
    }

    for (const auto& elt : mesh_file.elements) {
        if (elt.first < 0) {
            throw std::logic_error("ERROR in mesh_metadata.cpp: Element ID is negative; not supported\n");
        }

        uint ID = elt.first;
        _elements.insert({ID, ElementMetaData(3)});

        for (uint i = 0; i < 3; ++i) {
            _elements.at(ID).node_ids[i] = elt.second[i + 1];
        }
    }

    // make all edges
    using eltID_faceID = std::pair<uint, uint>;
    std::unordered_map<std::uint64_t, std::pair<eltID_faceID, eltID_faceID>> edge_dictionary;

    for (const auto& elt : mesh_file.elements) {
        std::vector<int> node{elt.second[1], elt.second[2], elt.second[3]};

        for (uint k = 0; k < 3; ++k) {
            std::uint64_t curr_key = (std::uint64_t)(std::min(node[(k + 1) % 3], node[(k + 2) % 3])) << 32 |
                                     std::max(node[(k + 1) % 3], node[(k + 2) % 3]);

            if (edge_dictionary.count(curr_key)) {  // if already one element has the edge.
                edge_dictionary.at(curr_key).second = std::make_pair(elt.first, k);
            } else {
                std::pair<eltID_faceID, eltID_faceID> edge_info{{elt.first, k}, {DEFAULT_ID, 0}};
                edge_dictionary.insert({curr_key, edge_info});
            }
        }
    }

    for (const auto& edge : edge_dictionary) {
        // check if there are two elements associated with this edge
        if (edge.second.second.first != DEFAULT_ID) {
            uint eltA_id = edge.second.first.first;
            uint eltB_id = edge.second.second.first;
            uint faceA_id = edge.second.first.second;
            uint faceB_id = edge.second.second.second;

            _elements.at(eltA_id).neighbor_ID[faceA_id] = eltB_id;
            _elements.at(eltA_id).boundary_type[faceA_id] = SWE::BoundaryConditions::internal;

            _elements.at(eltB_id).neighbor_ID[faceB_id] = eltA_id;
            _elements.at(eltB_id).boundary_type[faceB_id] = SWE::BoundaryConditions::internal;
        } else {
            // treat boundary conditions
            uint elt_id = edge.second.first.first;
            uint face_id = edge.second.first.second;

            std::array<int, 2> nodes{mesh_file.elements.at(elt_id)[(face_id + 1) % 3 + 1],
                                     mesh_file.elements.at(elt_id)[(face_id + 2) % 3 + 1]};

            _elements.at(elt_id).boundary_type[face_id] = mesh_file.get_ibtype(nodes);
        }
    }
}

std::vector<Point<2>> MeshMetaData::GetNodalCoordinates(uint elt_id) const {
    const std::vector<uint>& node_ids = _elements.at(elt_id).node_ids;
    std::vector<Point<2>> nodal_coordinates(node_ids.size());
    for (uint indx = 0; indx < node_ids.size(); ++indx) {
        nodal_coordinates[indx] = _nodes.at(node_ids[indx]).coordinates;
    }
    return nodal_coordinates;
}

std::vector<double> MeshMetaData::GetBathymetry(uint elt_id) const {
    const std::vector<uint>& node_ids = _elements.at(elt_id).node_ids;
    std::vector<double> nodal_bathymetry(node_ids.size());
    for (uint indx = 0; indx < node_ids.size(); ++indx) {
        nodal_bathymetry[indx] = _nodes.at(node_ids[indx]).bathymetry;
    }
    return nodal_bathymetry;
}