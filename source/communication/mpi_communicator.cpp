#include "mpi_communicator.hpp"
#include "preprocessor/mesh_metadata.hpp"

MPICommunicator::MPICommunicator(const std::string& neighborhood_data_file,
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

        MPIRankBoundary rank_boundary;

        std::string my_location;
        std::string neighbor_location;

        if (locality_A == locality_id && submesh_A == submesh_id) {
            rank_boundary.send_rank = locality_B;
            rank_boundary.receive_rank = locality_A;

            rank_boundary.send_tag = (int)((unsigned short)submesh_A << 16 | (unsigned short)submesh_B);
            rank_boundary.receive_tag = (int)((unsigned short)submesh_B << 16 | (unsigned short)submesh_A);
        } else {
            rank_boundary.send_rank = locality_A;
            rank_boundary.receive_rank = locality_B;

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

    this->send_requests.resize(this->GetRankBoundaryNumber());
    this->receieve_requests.resize(this->GetRankBoundaryNumber());
    
    this->send_statuses.resize(this->GetRankBoundaryNumber());
    this->receieve_statuses.resize(this->GetRankBoundaryNumber());
}

void MPICommunicator::SendAll(const uint timestamp) {
    for (uint rank_boundary_id = 0; rank_boundary_id < this->GetRankBoundaryNumber(); rank_boundary_id++) {
        this->send_requests[rank_boundary_id] = this->rank_boundaries[rank_boundary_id].send(timestamp);
    }
}

void MPICommunicator::ReceiveAll(const uint timestamp) {
    for (uint rank_boundary_id = 0; rank_boundary_id < this->GetRankBoundaryNumber(); rank_boundary_id++) {
        this->receive_requests[rank_boundary_id] = this->rank_boundaries[rank_boundary_id].receive(timestamp);
    }
}

void MPICommunicator::WaitAllSends(const uint timestamp) {
    MPI_Waitall(this->GetRankBoundaryNumber(), &this->send_requests.front(), &this->send_statuses.front());
}

void MPICommunicator::WaitAllReceives(const uint timestamp) {
    MPI_Waitall(this->GetRankBoundaryNumber(), &this->receive_requests.front(), &this->receive_statuses.front());
}        