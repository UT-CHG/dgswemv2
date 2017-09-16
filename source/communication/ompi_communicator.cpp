#include "ompi_communicator.hpp"

OMPICommunicator::OMPICommunicator(const std::string& neighborhood_data_file,
                                   const uint locality_id,
                                   const uint submesh_id) {
    std::ifstream file(neighborhood_data_file);

    if (!file) {
        throw std::logic_error("Error: Unable to find distributed boundary file : " + neighborhood_data_file + '\n');
    }

    std::string line;

    uint locality_A, locality_B, submesh_A, submesh_B, n_dboubdaries;

    while (std::getline(file, line)) {
        std::stringstream neighborhood_data(line);

        neighborhood_data >> locality_A >> submesh_A >> locality_B >> submesh_B >> n_dboubdaries;

        OMPIRankBoundary rank_boundary;

        if (locality_A == locality_id && submesh_A == submesh_id) {
            rank_boundary.send_rank = locality_B;
            rank_boundary.receive_rank = locality_B;

            rank_boundary.send_tag = (int)((unsigned short)submesh_A << 16 | (unsigned short)submesh_B);
            rank_boundary.receive_tag = (int)((unsigned short)submesh_B << 16 | (unsigned short)submesh_A);
        } else {
            rank_boundary.send_rank = locality_A;
            rank_boundary.receive_rank = locality_A;

            rank_boundary.send_tag = (int)((unsigned short)submesh_B << 16 | (unsigned short)submesh_A);
            rank_boundary.receive_tag = (int)((unsigned short)submesh_A << 16 | (unsigned short)submesh_B);
        }

        rank_boundary.elements.reserve(n_dboubdaries);
        rank_boundary.bound_ids.reserve(n_dboubdaries);
        rank_boundary.p.reserve(n_dboubdaries);

        for (uint db = 0; db < n_dboubdaries; ++db) {
            DistributedBoundaryMetaData dboundary_meta_data;

            file >> dboundary_meta_data;

            file.ignore(1000, '\n');

            if (locality_A == locality_id && submesh_A == submesh_id) {
                rank_boundary.elements.push_back(dboundary_meta_data.elements.first);
                rank_boundary.bound_ids.push_back(dboundary_meta_data.bound_ids.first);
                rank_boundary.p.push_back(dboundary_meta_data.p);
            } else {
                rank_boundary.elements.push_back(dboundary_meta_data.elements.second);
                rank_boundary.bound_ids.push_back(dboundary_meta_data.bound_ids.second);
                rank_boundary.p.push_back(dboundary_meta_data.p);
            }
        }

        this->rank_boundaries.push_back(std::move(rank_boundary));
    }

    this->send_requests.resize(this->rank_boundaries.size());
    this->receive_requests.resize(this->rank_boundaries.size());
}

void OMPICommunicator::InitializeCommunication() {
    for (uint rank_boundary_id = 0; rank_boundary_id < this->rank_boundaries.size(); rank_boundary_id++) {
        OMPIRankBoundary& rank_boundary = this->rank_boundaries[rank_boundary_id];

        MPI_Send_init(&rank_boundary.send_buffer.front(),
                      rank_boundary.send_buffer.size(),
                      MPI_DOUBLE,
                      rank_boundary.send_rank,
                      rank_boundary.send_tag,
                      MPI_COMM_WORLD,
                      &this->send_requests.at(rank_boundary_id));

        MPI_Recv_init(&rank_boundary.receive_buffer.front(),
                      rank_boundary.receive_buffer.size(),
                      MPI_DOUBLE,
                      rank_boundary.receive_rank,
                      rank_boundary.receive_tag,
                      MPI_COMM_WORLD,
                      &this->receive_requests.at(rank_boundary_id));
    }
}
