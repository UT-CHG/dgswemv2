#include "mesh_metadata.hpp"

MeshMetaData::MeshMetaData(const AdcircFormat& mesh_file) {
    _mesh_name = mesh_file.name;

    for (const auto& nod : mesh_file.nodes) {
        if (nod.first < 0) {
            throw std::logic_error("ERROR in mesh_metadata.cpp: Node ID is negative; not supported\n");
        }

        uint ID = nod.first;
        _nodes.insert({ID, {nod.second[0], nod.second[1], nod.second[2]}});
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
            std::uint64_t curr_key = static_cast<std::uint64_t>(std::min(node[(k + 1) % 3], node[(k + 2) % 3])) << 32 |
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

            _elements.at(elt_id).neighbor_ID[face_id] = DEFAULT_ID;
            _elements.at(elt_id).boundary_type[face_id] = mesh_file.get_ibtype(nodes);
        }
    }
}

MeshMetaData::MeshMetaData(const std::string& file) {
    std::ifstream ifs(file);

    if (!ifs) {
        std::string err_msg = "Fatal Error: Mesh named " + file + " not found\n";
        throw std::logic_error(err_msg);
    }

    std::getline(ifs, _mesh_name);

    uint num_elements;
    ifs >> num_elements;
    ifs.ignore(1000, '\n');
    uint elt_id;
    for (uint e = 0; e < num_elements; ++e) {
        ifs >> elt_id;
        ifs >> _elements[elt_id];
        ifs.ignore(1000, '\n');
    }

    uint num_nodes;
    ifs >> num_nodes;
    ifs.ignore(1000, '\n');
    uint node_id;
    for (uint n = 0; n < num_nodes; ++n) {
        ifs >> node_id;
        ifs >> _nodes[node_id];
        ifs.ignore(1000, '\n');
    }

    ifs.close();
}

void MeshMetaData::WriteTo(const std::string& file) {
    std::ofstream ofs;
    ofs.open(file);

    ofs << _mesh_name << '\n';
    ofs << _elements.size() << " = number of elements\n";
    for (const auto& elt : _elements) {
        ofs << elt.first << " " << elt.second << '\n';
    }

    ofs << _nodes.size() << " = number of nodes\n";
    for (const auto& nod : _nodes) {
        ofs << nod.first << " " << nod.second << '\n';
    }

    ofs.close();
}

std::vector<Point<3>> MeshMetaData::GetNodalCoordinates(uint elt_id) const {
    const std::vector<uint>& node_ids = _elements.at(elt_id).node_ids;

    std::vector<Point<3>> nodal_coordinates(node_ids.size());

    for (uint indx = 0; indx < node_ids.size(); ++indx) {
        nodal_coordinates[indx] = _nodes.at(node_ids[indx]).coordinates;
    }

    return nodal_coordinates;
}
