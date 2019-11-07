#include "ompi_communicator.hpp"

OMPICommunicator::OMPICommunicator(const DistributedBoundaryMetaData& db_data) {
    for (auto& rb_meta_data : db_data.rank_boundary_data) {
        OMPIRankBoundary rank_boundary;

        rank_boundary.db_data = rb_meta_data;

        rank_boundary.send_rank    = rb_meta_data.locality_ex;
        rank_boundary.receive_rank = rb_meta_data.locality_ex;

        rank_boundary.send_tag =
            (int)((unsigned short)rb_meta_data.submesh_in << 8 | (unsigned short)rb_meta_data.submesh_ex);
        rank_boundary.receive_tag =
            (int)((unsigned short)rb_meta_data.submesh_ex << 8 | (unsigned short)rb_meta_data.submesh_in);

        this->rank_boundaries.push_back(std::move(rank_boundary));
    }
}

void OMPICommunicator::InitializeCommunication() {
    uint ncomm   = this->rank_boundaries.begin()->send_buffer.size();
    uint nrbound = this->rank_boundaries.size();

    this->send_requests.resize(ncomm);
    this->receive_requests.resize(ncomm);

    for (uint comm = 0; comm < ncomm; ++comm) {
        this->send_requests[comm].resize(nrbound);
        this->receive_requests[comm].resize(nrbound);
    }

    for (uint rank_boundary_id = 0; rank_boundary_id < this->rank_boundaries.size(); ++rank_boundary_id) {
        OMPIRankBoundary& rank_boundary = this->rank_boundaries[rank_boundary_id];

        uint ncomm = rank_boundary.send_buffer.size();

        for (uint comm = 0; comm < ncomm; ++comm) {
            MPI_Request& send_request    = this->send_requests[comm][rank_boundary_id];
            MPI_Request& receive_request = this->receive_requests[comm][rank_boundary_id];

            MPI_Send_init(&rank_boundary.send_buffer[comm].front(),
                          rank_boundary.send_buffer[comm].size(),
                          MPI_DOUBLE,
                          rank_boundary.send_rank,
                          rank_boundary.send_tag,
                          MPI_COMM_WORLD,
                          &send_request);

            MPI_Recv_init(&rank_boundary.receive_buffer[comm].front(),
                          rank_boundary.receive_buffer[comm].size(),
                          MPI_DOUBLE,
                          rank_boundary.receive_rank,
                          rank_boundary.receive_tag,
                          MPI_COMM_WORLD,
                          &receive_request);
        }
    }
}

void OMPICommunicator::SendAll(const uint comm_type, const uint timestamp) {
    MPI_Startall(this->send_requests[comm_type].size(), &this->send_requests[comm_type].front());
}

void OMPICommunicator::ReceiveAll(const uint comm_type, const uint timestamp) {
    MPI_Startall(this->receive_requests[comm_type].size(), &this->receive_requests[comm_type].front());
}

void OMPICommunicator::WaitAllSends(const uint comm_type, const uint timestamp) {
    MPI_Waitall(this->send_requests[comm_type].size(), &this->send_requests[comm_type].front(), MPI_STATUSES_IGNORE);
}

void OMPICommunicator::WaitAllReceives(const uint comm_type, const uint timestamp) {
    MPI_Waitall(
        this->receive_requests[comm_type].size(), &this->receive_requests[comm_type].front(), MPI_STATUSES_IGNORE);
}