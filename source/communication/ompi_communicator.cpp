#include "ompi_communicator.hpp"

OMPICommunicator::OMPICommunicator(const DistributedBoundaryMetaData& db_data) {
    for (auto& rb_meta_data : db_data.rank_boundary_data) {
        OMPIRankBoundary rank_boundary;

        rank_boundary.db_data = rb_meta_data;

        rank_boundary.send_rank    = rb_meta_data.locality_ex;
        rank_boundary.receive_rank = rb_meta_data.locality_ex;

        rank_boundary.send_tag =
            (int)((unsigned short)rb_meta_data.submesh_in << 16 | (unsigned short)rb_meta_data.submesh_ex);
        rank_boundary.receive_tag =
            (int)((unsigned short)rb_meta_data.submesh_ex << 16 | (unsigned short)rb_meta_data.submesh_in);

        this->rank_boundaries.push_back(std::move(rank_boundary));
    }

    this->send_preproc_requests.resize(this->rank_boundaries.size());
    this->receive_preproc_requests.resize(this->rank_boundaries.size());

    this->send_requests.resize(this->rank_boundaries.size());
    this->receive_requests.resize(this->rank_boundaries.size());

    this->send_postproc_requests.resize(this->rank_boundaries.size());
    this->receive_postproc_requests.resize(this->rank_boundaries.size());
}

void OMPICommunicator::InitializeCommunication() {
    for (uint rank_boundary_id = 0; rank_boundary_id < this->rank_boundaries.size(); rank_boundary_id++) {
        OMPIRankBoundary& rank_boundary = this->rank_boundaries[rank_boundary_id];

        MPI_Send_init(&rank_boundary.send_preproc_buffer.front(),
                      rank_boundary.send_preproc_buffer.size(),
                      MPI_DOUBLE,
                      rank_boundary.send_rank,
                      rank_boundary.send_tag,
                      MPI_COMM_WORLD,
                      &this->send_preproc_requests.at(rank_boundary_id));

        MPI_Recv_init(&rank_boundary.receive_preproc_buffer.front(),
                      rank_boundary.receive_preproc_buffer.size(),
                      MPI_DOUBLE,
                      rank_boundary.receive_rank,
                      rank_boundary.receive_tag,
                      MPI_COMM_WORLD,
                      &this->receive_preproc_requests.at(rank_boundary_id));

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

        MPI_Send_init(&rank_boundary.send_postproc_buffer.front(),
                      rank_boundary.send_postproc_buffer.size(),
                      MPI_DOUBLE,
                      rank_boundary.send_rank,
                      rank_boundary.send_tag,
                      MPI_COMM_WORLD,
                      &this->send_postproc_requests.at(rank_boundary_id));

        MPI_Recv_init(&rank_boundary.receive_postproc_buffer.front(),
                      rank_boundary.receive_postproc_buffer.size(),
                      MPI_DOUBLE,
                      rank_boundary.receive_rank,
                      rank_boundary.receive_tag,
                      MPI_COMM_WORLD,
                      &this->receive_postproc_requests.at(rank_boundary_id));
    }
}

void OMPICommunicator::SendPreprocAll(const uint timestamp) {
    MPI_Startall(this->send_preproc_requests.size(), &this->send_preproc_requests.front());
}

void OMPICommunicator::ReceivePreprocAll(const uint timestamp) {
    MPI_Startall(this->receive_preproc_requests.size(), &this->receive_preproc_requests.front());
}

void OMPICommunicator::WaitAllPreprocSends(const uint timestamp) {
    MPI_Waitall(this->send_preproc_requests.size(), &this->send_preproc_requests.front(), MPI_STATUSES_IGNORE);
}

void OMPICommunicator::WaitAllPreprocReceives(const uint timestamp) {
    MPI_Waitall(this->receive_preproc_requests.size(), &this->receive_preproc_requests.front(), MPI_STATUSES_IGNORE);
}

void OMPICommunicator::SendAll(const uint timestamp) {
    MPI_Startall(this->send_requests.size(), &this->send_requests.front());
}

void OMPICommunicator::ReceiveAll(const uint timestamp) {
    MPI_Startall(this->receive_requests.size(), &this->receive_requests.front());
}

void OMPICommunicator::WaitAllSends(const uint timestamp) {
    MPI_Waitall(this->send_requests.size(), &this->send_requests.front(), MPI_STATUSES_IGNORE);
}

void OMPICommunicator::WaitAllReceives(const uint timestamp) {
    MPI_Waitall(this->receive_requests.size(), &this->receive_requests.front(), MPI_STATUSES_IGNORE);
}

void OMPICommunicator::SendPostprocAll(const uint timestamp) {
    MPI_Startall(this->send_postproc_requests.size(), &this->send_postproc_requests.front());
}

void OMPICommunicator::ReceivePostprocAll(const uint timestamp) {
    MPI_Startall(this->receive_postproc_requests.size(), &this->receive_postproc_requests.front());
}

void OMPICommunicator::WaitAllPostprocSends(const uint timestamp) {
    MPI_Waitall(this->send_postproc_requests.size(), &this->send_postproc_requests.front(), MPI_STATUSES_IGNORE);
}

void OMPICommunicator::WaitAllPostprocReceives(const uint timestamp) {
    MPI_Waitall(this->receive_postproc_requests.size(), &this->receive_postproc_requests.front(), MPI_STATUSES_IGNORE);
}
