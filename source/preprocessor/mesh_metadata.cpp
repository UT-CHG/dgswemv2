#include "mesh_metadata.hpp"

MeshMetaData::MeshMetaData(const AdcircFormat& mesh_file) {
    this->mesh_name = mesh_file.name;

    for (const auto& nod : mesh_file.nodes) {
        if (nod.first < 0) {
            throw std::logic_error("Fatal Error: in ADCIRC mesh node ID is negative!\n");
        }

        uint ID = nod.first;
        this->nodes.insert({ID, {nod.second[0], nod.second[1], nod.second[2]}});
    }

    for (const auto& elt : mesh_file.elements) {
        if (elt.first < 0) {
            throw std::logic_error("Fatal Error: in ADCIRC mesh element ID is negative!\n");
        }

        uint ID = elt.first;
        this->elements.insert({ID, ElementMetaData(3)});

        for (uint i = 0; i < 3; ++i) {
            this->elements.at(ID).node_ID[i] = elt.second[i + 1];
        }
    }

    // make all edges
    using eltID_faceID = std::pair<uint, uint>;
    std::unordered_map<std::uint64_t, std::pair<eltID_faceID, eltID_faceID>> edge_dictionary;

    for (const auto& elt : mesh_file.elements) {
        std::vector<uint> node{elt.second[1], elt.second[2], elt.second[3]};

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

            this->elements.at(eltA_id).neighbor_ID[faceA_id]   = eltB_id;
            this->elements.at(eltA_id).boundary_type[faceA_id] = SWE::BoundaryTypes::internal;

            this->elements.at(eltB_id).neighbor_ID[faceB_id]   = eltA_id;
            this->elements.at(eltB_id).boundary_type[faceB_id] = SWE::BoundaryTypes::internal;
        } else {
            // treat boundary conditions
            uint elt_id  = edge.second.first.first;
            uint face_id = edge.second.first.second;

            std::array<uint, 2> nodes{mesh_file.elements.at(elt_id)[(face_id + 1) % 3 + 1],
                                      mesh_file.elements.at(elt_id)[(face_id + 2) % 3 + 1]};

            this->elements.at(elt_id).neighbor_ID[face_id]   = DEFAULT_ID;
            this->elements.at(elt_id).boundary_type[face_id] = mesh_file.get_ibtype(nodes);

            // find element id on the other side of internal barrier
            if (this->elements.at(elt_id).boundary_type[face_id] == SWE::BoundaryTypes::levee) {
                std::array<uint, 2> barrier_np = mesh_file.get_barrier_node_pair(nodes);
                std::uint64_t key = static_cast<std::uint64_t>(std::min(barrier_np[0], barrier_np[1])) << 32 |
                                    std::max(barrier_np[0], barrier_np[1]);

                this->elements.at(elt_id).neighbor_ID[face_id] = edge_dictionary.at(key).first.first;
            }
        }
    }
}

MeshMetaData::MeshMetaData(const std::string& mesh_file) {
    if (!Utilities::file_exists(mesh_file)) {
        throw std::logic_error("Fatal Error: meta mesh file " + mesh_file + " was not found!\n");
    }

    std::ifstream ifs(mesh_file);

    ifs >> this->mesh_name;
    ifs.ignore(1000, '\n');

    uint numelements;
    ifs >> numelements;
    ifs.ignore(1000, '\n');
    uint elt_id;
    for (uint e = 0; e < numelements; ++e) {
        ifs >> elt_id;
        ifs >> this->elements[elt_id];
        ifs.ignore(1000, '\n');
    }

    uint numnodes;
    ifs >> numnodes;
    ifs.ignore(1000, '\n');
    uint node_id;
    for (uint n = 0; n < numnodes; ++n) {
        ifs >> node_id;
        ifs >> this->nodes[node_id];
        ifs.ignore(1000, '\n');
    }

    ifs.close();
}

void MeshMetaData::write_to(const std::string& file) {
    std::ofstream ofs;
    ofs.open(file);

    ofs << this->mesh_name << '\n';
    ofs << this->elements.size() << " = number of elements\n";
    for (const auto& elt : this->elements) {
        ofs << elt.first << " " << elt.second << '\n';
    }

    ofs << this->nodes.size() << " = number of nodes\n";
    for (const auto& nod : this->nodes) {
        ofs << nod.first << " " << nod.second << '\n';
    }

    ofs.close();
}

std::vector<Point<3>> MeshMetaData::get_nodal_coordinates(uint elt_id) const {
    const std::vector<uint>& node_ID = this->elements.at(elt_id).node_ID;

    std::vector<Point<3>> nodal_coordinates(node_ID.size());

    for (uint indx = 0; indx < node_ID.size(); ++indx) {
        nodal_coordinates[indx] = this->nodes.at(node_ID[indx]).coordinates;
    }

    return nodal_coordinates;
}

DistributedBoundaryMetaData::DistributedBoundaryMetaData(const std::string& dbmd_file,
                                                         uint locality_id,
                                                         uint submesh_id) {
    if (!Utilities::file_exists(dbmd_file)) {
        throw std::logic_error("Fatal Error: distributed boundary data file " + dbmd_file + " was not found!\n");
    }

    std::ifstream file(dbmd_file);

    std::string line;

    uint locality_A, locality_B, submesh_A, submesh_B, n_dboubdaries;

    while (std::getline(file, line)) {
        std::stringstream neighborhood_data(line);

        neighborhood_data >> locality_A >> submesh_A >> locality_B >> submesh_B >> n_dboubdaries;

        RankBoundaryMetaData rank_boundary;

        if (locality_A == locality_id && submesh_A == submesh_id) {
            rank_boundary.locality_in = locality_A;
            rank_boundary.locality_ex = locality_B;

            rank_boundary.submesh_in = submesh_A;
            rank_boundary.submesh_ex = submesh_B;
        } else if (locality_B == locality_id && submesh_B == submesh_id) {
            rank_boundary.locality_in = locality_B;
            rank_boundary.locality_ex = locality_A;

            rank_boundary.submesh_in = submesh_B;
            rank_boundary.submesh_ex = submesh_A;
        } else {
            throw std::logic_error("Fatal Error: error in locality/submesh in distributed boundary file " + dbmd_file +
                                   "!\n");
        }

        rank_boundary.elements_in.reserve(n_dboubdaries);
        rank_boundary.elements_ex.reserve(n_dboubdaries);

        rank_boundary.bound_ids_in.reserve(n_dboubdaries);
        rank_boundary.bound_ids_ex.reserve(n_dboubdaries);

        rank_boundary.p.reserve(n_dboubdaries);

        for (uint db = 0; db < n_dboubdaries; ++db) {
            DBPairMetaData dboundary_meta_data;

            file >> dboundary_meta_data;

            file.ignore(1000, '\n');

            if (locality_A == locality_id && submesh_A == submesh_id) {
                rank_boundary.elements_in.push_back(dboundary_meta_data.elements.first);
                rank_boundary.elements_ex.push_back(dboundary_meta_data.elements.second);

                rank_boundary.bound_ids_in.push_back(dboundary_meta_data.bound_ids.first);
                rank_boundary.bound_ids_ex.push_back(dboundary_meta_data.bound_ids.second);

                rank_boundary.p.push_back(dboundary_meta_data.p);
            } else {
                rank_boundary.elements_in.push_back(dboundary_meta_data.elements.second);
                rank_boundary.elements_ex.push_back(dboundary_meta_data.elements.first);

                rank_boundary.bound_ids_in.push_back(dboundary_meta_data.bound_ids.second);
                rank_boundary.bound_ids_ex.push_back(dboundary_meta_data.bound_ids.first);

                rank_boundary.p.push_back(dboundary_meta_data.p);
            }
        }

        this->rank_boundary_data.push_back(std::move(rank_boundary));
    }
}
