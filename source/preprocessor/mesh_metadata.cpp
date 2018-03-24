#include "mesh_metadata.hpp"

MeshMetaData::MeshMetaData(const AdcircFormat& mesh_file) {
    mesh_name = mesh_file.name;

    for (const auto& nod : mesh_file.nodes) {
        if (nod.first < 0) {
            throw std::logic_error(
                "ERROR in mesh_metadata.cpp: Node ID is negative; not "
                "supported\n");
        }

        uint ID = nod.first;
        nodes.insert({ID, {nod.second[0], nod.second[1], nod.second[2]}});
    }

    for (const auto& elt : mesh_file.elements) {
        if (elt.first < 0) {
            throw std::logic_error(
                "ERROR in mesh_metadata.cpp: Element ID is negative; not "
                "supported\n");
        }

        uint ID = elt.first;
        elements.insert({ID, ElementMetaData(3)});

        for (uint i = 0; i < 3; ++i) {
            elements.at(ID).node_ID[i] = elt.second[i + 1];
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
            uint eltA_id  = edge.second.first.first;
            uint eltB_id  = edge.second.second.first;
            uint faceA_id = edge.second.first.second;
            uint faceB_id = edge.second.second.second;

            elements.at(eltA_id).neighbor_ID[faceA_id]   = eltB_id;
            elements.at(eltA_id).boundary_type[faceA_id] = INTERNAL;

            elements.at(eltB_id).neighbor_ID[faceB_id]   = eltA_id;
            elements.at(eltB_id).boundary_type[faceB_id] = INTERNAL;
        } else {
            // treat boundary conditions
            uint elt_id  = edge.second.first.first;
            uint face_id = edge.second.first.second;

            std::array<int, 2> nodes{mesh_file.elements.at(elt_id)[(face_id + 1) % 3 + 1],
                                     mesh_file.elements.at(elt_id)[(face_id + 2) % 3 + 1]};

            elements.at(elt_id).neighbor_ID[face_id]   = DEFAULT_ID;
            elements.at(elt_id).boundary_type[face_id] = mesh_file.get_ibtype(nodes);
        }
    }
}

MeshMetaData::MeshMetaData(const std::string& file) {
    std::ifstream ifs(file);

    if (!ifs) {
        std::string err_msg = "Fatal Error: Mesh named " + file + " not found\n";
        throw std::logic_error(err_msg);
    }

    std::getline(ifs, mesh_name);

    uint numelements;
    ifs >> numelements;
    ifs.ignore(1000, '\n');
    uint elt_id;
    for (uint e = 0; e < numelements; ++e) {
        ifs >> elt_id;
        ifs >> elements[elt_id];
        ifs.ignore(1000, '\n');
    }

    uint numnodes;
    ifs >> numnodes;
    ifs.ignore(1000, '\n');
    uint node_id;
    for (uint n = 0; n < numnodes; ++n) {
        ifs >> node_id;
        ifs >> nodes[node_id];
        ifs.ignore(1000, '\n');
    }

    ifs.close();
}

void MeshMetaData::write_to(const std::string& file) {
    std::ofstream ofs;
    ofs.open(file);

    ofs << mesh_name << '\n';
    ofs << elements.size() << " = number of elements\n";
    for (const auto& elt : elements) {
        ofs << elt.first << " " << elt.second << '\n';
    }

    ofs << nodes.size() << " = number of nodes\n";
    for (const auto& nod : nodes) {
        ofs << nod.first << " " << nod.second << '\n';
    }

    ofs.close();
}

std::vector<Point<3>> MeshMetaData::get_nodal_coordinates(uint elt_id) const {
    const std::vector<uint>& node_ID = elements.at(elt_id).node_ID;

    std::vector<Point<3>> nodal_coordinates(node_ID.size());

    for (uint indx = 0; indx < node_ID.size(); ++indx) {
        nodal_coordinates[indx] = nodes.at(node_ID[indx]).coordinates;
    }

    return nodal_coordinates;
}
