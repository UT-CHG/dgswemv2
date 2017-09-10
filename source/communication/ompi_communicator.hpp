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
};

class OMPICommunicator {
  private:
    std::vector<OMPIRankBoundary> rank_boundaries;
    
    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> receive_requests;
    
  public:
    OMPICommunicator() = default;
    OMPICommunicator(const std::string& neighborhood_data_file, const uint locality_id, const uint submesh_id);

    void InitializeCommunication();

    uint GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    OMPIRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }
    
    void SendAll(const uint timestamp) {
        MPI_Startall(this->send_requests.size(), &this->send_requests.front(), MPI_STATUSES_IGNORE)
    }
  
    void ReceiveAll(const uint timestamp) {
       MPI_Startall(this->receive_requests.size(), &this->receive_requests.front(), MPI_STATUSES_IGNORE)
    }
  
    void WaitAllSends(const uint timestamp) {
      MPI_Waitall(this->send_requests.size(), &this->send_requests.front(), MPI_STATUSES_IGNORE)
    }

    void WaitAllReceives(const uint timestamp) {
      MPI_Waitall(this->receive_requests.size(), &this->receive_requests.front(), MPI_STATUSES_IGNORE)
    }

  public:
    using RankBoundaryType = OMPIRankBoundary;
};

#endif
