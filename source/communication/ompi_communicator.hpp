#ifndef OMPI_COMMUNICATOR_HPP
#define OMPI_COMMUNICATOR_HPP

#include <mpi.h>

#include "general_definitions.hpp"
#include "preprocessor/mesh_metadata.hpp"

struct OMPIRankBoundary {
    RankBoundaryMetaData db_data;

    int send_rank;
    int receive_rank;

    int send_tag;
    int receive_tag;

    std::vector<std::vector<double>> send_buffer;
    std::vector<std::vector<double>> receive_buffer;
};

class OMPICommunicator {
  private:
    std::vector<OMPIRankBoundary> rank_boundaries;

    std::vector<std::vector<MPI_Request>> send_requests;
    std::vector<std::vector<MPI_Request>> receive_requests;

  public:
    OMPICommunicator() = default;
    OMPICommunicator(const DistributedBoundaryMetaData& db_data);

    void InitializeCommunication();

    uint GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    OMPIRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }

    void SendAll(const uint comm_type, const uint timestamp);
    void ReceiveAll(const uint comm_type, const uint timestamp);
    void WaitAllSends(const uint comm_type, const uint timestamp);
    void WaitAllReceives(const uint comm_type, const uint timestamp);

  public:
    using RankBoundaryType = OMPIRankBoundary;
};

#endif
