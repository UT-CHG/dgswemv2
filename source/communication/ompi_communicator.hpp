#ifndef OMPI_COMMUNICATOR_HPP
#define OMPI_COMMUNICATOR_HPP

#include <mpi.h>

#include "../general_definitions.hpp"
#include "../preprocessor/mesh_metadata.hpp"

struct OMPIRankBoundary {
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
    OMPICommunicator(const DistributedBoundaryMetaData& db_data);

    void InitializeCommunication();

    uint GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    OMPIRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }

    void SendPreprocAll(const uint timestamp);
    void ReceivePreprocAll(const uint timestamp);
    void WaitAllPreprocSends(const uint timestamp);
    void WaitAllPreprocReceives(const uint timestamp);

    void SendAll(const uint timestamp);
    void ReceiveAll(const uint timestamp);
    void WaitAllSends(const uint timestamp);
    void WaitAllReceives(const uint timestamp);

    void SendPostprocAll(const uint timestamp);
    void ReceivePostprocAll(const uint timestamp);
    void WaitAllPostprocSends(const uint timestamp);
    void WaitAllPostprocReceives(const uint timestamp);

  public:
    using RankBoundaryType = OMPIRankBoundary;
};

#endif
