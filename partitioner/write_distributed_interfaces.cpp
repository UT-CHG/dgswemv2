#include "preprocessor/mesh_metadata.hpp"
#include "preprocessor/input_parameters.hpp"
#include "util.hpp"

#include <deque>
#include <unordered_set>

void write_distributed_edge_metadata(const std::string& file_name,
                                     const InputParameters& input,
                                     const MeshMetaData& mesh_meta,
                                     const std::vector<std::vector<MeshMetaData>>& submeshes) {
    std::size_t num_loc = submeshes.size();

    std::unordered_set<std::pair<uint, uint>> faces;
    // assemble faces shared across elements
    for (auto& elt : mesh_meta._elements) {
        for (auto& neigh : elt.second.neighbor_ID) {
            if (neigh != DEFAULT_ID) {
                std::pair<uint, uint> edge_name{std::min(elt.first, neigh), std::max(elt.first, neigh)};
                faces.insert(edge_name);
            }
        }
    }
    std::unordered_map<uint, uint> elt2partition;
    for (uint loc_id = 0; loc_id < submeshes.size(); ++loc_id) {
        for (uint sbmsh_id = 0; sbmsh_id < submeshes[loc_id].size(); ++sbmsh_id) {
            for (auto& elt : submeshes[loc_id][sbmsh_id]._elements) {
                elt2partition.insert(std::make_pair(elt.first, loc_id + num_loc * sbmsh_id));
                for (auto& neigh : elt.second.neighbor_ID) {
                    if (neigh != DEFAULT_ID) {
                        std::pair<uint, uint> edge_name{std::min(elt.first, neigh), std::max(elt.first, neigh)};
                        faces.erase(edge_name);
                    }
                }
            }
        }
    }
    // Assemble lists of distributed interface meta data
    std::unordered_map<std::pair<uint, uint>, std::deque<DistributedInterfaceMetaData>> shared_faces;
    for (auto& f : faces) {
        uint eltA = f.first;
        uint rnkA = elt2partition.at(eltA);

        uint eltB = f.second;
        uint rnkB = elt2partition.at(eltB);

        uint face_id_A{DEFAULT_ID};
        const ElementMetaData& eltA_meta = mesh_meta._elements.at(eltA);
        for (uint fid = 0; fid < eltA_meta.neighbor_ID.size(); ++fid) {
            if (eltB == eltA_meta.neighbor_ID[fid]) {
                face_id_A = fid;
            }
        }
        assert(face_id_A != DEFAULT_ID);

        uint face_id_B{DEFAULT_ID};
        const ElementMetaData& eltB_meta = mesh_meta._elements.at(eltB);
        for (uint fid = 0; fid < eltB_meta.neighbor_ID.size(); ++fid) {
            if (eltA == eltB_meta.neighbor_ID[fid]) {
                face_id_B = fid;
            }
        }
        assert(face_id_B != DEFAULT_ID);

        std::pair<uint, uint> rnk_pair{std::min(rnkA, rnkB), std::max(rnkA, rnkB)};

        DistributedInterfaceMetaData dist_int;
        dist_int.polynomial_order = input.polynomial_order;

        if (eltA > eltB) {
            dist_int.elements = {eltB, eltA};
            dist_int.face_id = {face_id_B, face_id_A};
        } else {
            dist_int.elements = {eltA, eltB};
            dist_int.face_id = {face_id_A, face_id_B};
        }

        shared_faces[rnk_pair].push_back(std::move(dist_int));
    }

    for (uint loc_id = 0; loc_id < submeshes.size(); ++loc_id) {
        for (uint sbmsh_id = 0; sbmsh_id < submeshes[loc_id].size(); ++sbmsh_id) {

            std::string distributed_meta_filename = file_name;
            distributed_meta_filename.erase(distributed_meta_filename.size() - 3);
            std::ofstream file(distributed_meta_filename + '_' + std::to_string(loc_id) + '_' +
                               std::to_string(sbmsh_id) + ".14_meta");
            for (auto& sf : shared_faces) {
                uint locA = sf.first.first % num_loc;
                uint sbmshA = sf.first.first / num_loc;

                uint locB = sf.first.second % num_loc;
                uint sbmshB = sf.first.second / num_loc;
                if ((locA == loc_id && sbmshA == sbmsh_id) || (locB == loc_id && sbmshB == sbmsh_id)) {
                    file << locA << " " << sbmshA << " " << locB << " " << sbmshB << " " << sf.second.size() << '\n';
                    for (auto& dist_int : sf.second) {
                        file << dist_int << '\n';
                    }
                }
            }
        }
    }
}