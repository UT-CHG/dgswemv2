#ifndef HPX_COMMUNICATOR_HPP
#define HPX_COMMUNICATOR_HPP

#include <hpx/hpx.hpp>
#include <vector>

using array_double = std::vector<double>;
HPX_REGISTER_CHANNEL_DECLARATION(array_double);

struct HPXRankBoundary {
    std::vector<uint> elements;
    std::vector<uint> bound_ids;
    std::vector<uint> p;

    hpx::lcos::channel<array_double> outgoing;
    hpx::lcos::channel<array_double> incoming;

    std::vector<double> send_buffer;
    std::vector<double> receive_buffer;

    void send(uint timestamp) { outgoing.set(send_buffer, timestamp); }

    hpx::future<void> receive(uint timestamp) {
        return incoming.get(timestamp)
            .then([this](hpx::future<array_double> msg_future) { this->receive_buffer = msg_future.get(); });
    }
};

class HPXCommunicator {
  private:
    std::vector<HPXRankBoundary> rank_boundaries;

  public:
    HPXCommunicator() = default;
    HPXCommunicator(const std::string&, uint, uint);

    void SendAll(uint);
    hpx::future<void> ReceiveAll(uint);

    uint GetRankBoundaryNumber() { return rank_boundaries.size(); }
    HPXRankBoundary& GetRankBoundary(uint rank_boundary_id) { return rank_boundaries.at(rank_boundary_id); }

  public:
    using RankBoundaryType = HPXRankBoundary;
};

#endif
