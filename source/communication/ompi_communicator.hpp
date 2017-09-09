#ifndef OMPI_COMMUNICATOR_HPP
#define OMPI_COMMUNICATOR_HPP

#include <mpi.h>
#include "../general_definitions.hpp"

struct OMPIRankBoundary {
    std::vector<uint> elements;
    std::vector<uint> bound_ids;
    std::vector<uint> p;

    int send_rank;
    int receive_rank;

    int send_tag;
    int receive_tag;

    std::vector<double> send_buffer;
    std::vector<double> receive_buffer;

    MPI_Request send(uint timestamp) {
        MPI_Request req;

        MPI_Isend(&this->send_buffer.front(),
                  this->send_buffer.size(),
                  MPI_DOUBLE,
                  this->send_rank,
                  this->send_tag,
                  MPI_COMM_WORLD,
                  &req);

        return req;
    }

    MPI_Request receive(uint timestamp) {
        MPI_Request req;

        MPI_Irecv(&this->receive_buffer.front(),
                  this->receive_buffer.size(),
                  MPI_DOUBLE,
                  this->receive_rank,
                  this->receive_tag,
                  MPI_COMM_WORLD,
                  &req);

        return req;
    }
};

class OMPICommunicator {
  private:
    std::vector<OMPIRankBoundary> rank_boundaries;

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> receive_requests;

    std::vector<MPI_Status> send_statuses;
    std::vector<MPI_Status> receive_statuses;

  public:
    OMPICommunicator() = default;
    OMPICommunicator(const std::string& neighborhood_data_file, const uint locality_id, const uint submesh_id);

    uint GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    OMPIRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }

    void SendAll(const uint timestamp);
    void ReceiveAll(const uint timestamp);

    void WaitAllSends(const uint timestamp);
    void WaitAllReceives(const uint timestamp);

  public:
    using RankBoundaryType = OMPIRankBoundary;
};

#endif
