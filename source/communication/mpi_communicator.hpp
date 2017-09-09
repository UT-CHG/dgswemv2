#ifndef MPI_COMMUNICATOR_HPP
#define MPI_COMMUNICATOR_HPP

struct MPIRankBoundary {
    std::vector<uint> elements;
    std::vector<uint> bound_ids;
    std::vector<uint> p;

    std::vector<double> send_buffer;
    std::vector<double> receive_buffer;

    void send(uint timestamp) { /*outgoing.set(send_buffer, timestamp);*/ }

    void receive(uint timestamp) {
        /*return incoming.get(timestamp)
            .then([this](hpx::future<array_double> msg_future) { this->receive_buffer = msg_future.get(); });*/
    }
};

class MPICommunicator {
  private:
    std::vector<MPIRankBoundary> rank_boundaries;

  public:
    MPICommunicator() = default;
    MPICommunicator(const std::string& neighborhood_data_file, const uint locality_id, const uint submesh_id);

    void SendAll(const uint timestamp);
    void ReceiveAll(const uint timestamp);

    uint GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    MPIRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }

  public:
    using RankBoundaryType = MPIRankBoundary;
};

#endif
