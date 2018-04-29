#ifndef OMPI_COMMUNICATOR_HPP
#define OMPI_COMMUNICATOR_HPP

#include <mpi.h>

#include "../general_definitions.hpp"
#include "../preprocessor/mesh_metadata.hpp"

struct OMPIRankBoundary {
    std::vector<uint> elements_in;
    std::vector<uint> elements_ex;

    std::vector<uint> bound_ids_in;
    std::vector<uint> bound_ids_ex;

    std::vector<uint> p;

    int send_rank;
    int receive_rank;

    int send_tag;
    int receive_tag;

    std::vector<double> send_preproc_buffer;
    std::vector<double> receive_preproc_buffer;

    std::vector<double> send_buffer;
    std::vector<double> receive_buffer;

    std::vector<double> send_postproc_buffer;
    std::vector<double> receive_postproc_buffer;
};

class OMPICommunicator {
  private:
    std::vector<OMPIRankBoundary> rank_boundaries;

    std::vector<MPI_Request> send_preproc_requests;
    std::vector<MPI_Request> receive_preproc_requests;

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> receive_requests;

    std::vector<MPI_Request> send_postproc_requests;
    std::vector<MPI_Request> receive_postproc_requests;

  public:
    OMPICommunicator() = default;
    OMPICommunicator(const std::string& neighborhood_data_file, const uint locality_id, const uint submesh_id);

    void InitializeCommunication();

    uint              GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    OMPIRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }

    void SendPreprocAll(const uint timestamp) {
        MPI_Startall(this->send_preproc_requests.size(), &this->send_preproc_requests.front());
    }

    void ReceivePreprocAll(const uint timestamp) {
        MPI_Startall(this->receive_preproc_requests.size(), &this->receive_preproc_requests.front());
    }

    void WaitAllPreprocSends(const uint timestamp) {
        MPI_Waitall(this->send_preproc_requests.size(), &this->send_preproc_requests.front(), MPI_STATUSES_IGNORE);
    }

    void WaitAllPreprocReceives(const uint timestamp) {
        MPI_Waitall(
            this->receive_preproc_requests.size(), &this->receive_preproc_requests.front(), MPI_STATUSES_IGNORE);
    }

    void SendAll(const uint timestamp) { MPI_Startall(this->send_requests.size(), &this->send_requests.front()); }

    void ReceiveAll(const uint timestamp) {
        MPI_Startall(this->receive_requests.size(), &this->receive_requests.front());
    }

    void WaitAllSends(const uint timestamp) {
        MPI_Waitall(this->send_requests.size(), &this->send_requests.front(), MPI_STATUSES_IGNORE);
    }

    void WaitAllReceives(const uint timestamp) {
        MPI_Waitall(this->receive_requests.size(), &this->receive_requests.front(), MPI_STATUSES_IGNORE);
    }

    void SendPostprocAll(const uint timestamp) {
        MPI_Startall(this->send_postproc_requests.size(), &this->send_postproc_requests.front());
    }

    void ReceivePostprocAll(const uint timestamp) {
        MPI_Startall(this->receive_postproc_requests.size(), &this->receive_postproc_requests.front());
    }

    void WaitAllPostprocSends(const uint timestamp) {
        MPI_Waitall(this->send_postproc_requests.size(), &this->send_postproc_requests.front(), MPI_STATUSES_IGNORE);
    }

    void WaitAllPostprocReceives(const uint timestamp) {
        MPI_Waitall(
            this->receive_postproc_requests.size(), &this->receive_postproc_requests.front(), MPI_STATUSES_IGNORE);
    }

  public:
    using RankBoundaryType = OMPIRankBoundary;
};

#endif
