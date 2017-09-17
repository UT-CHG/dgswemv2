#ifndef HPX_COMMUNICATOR_HPP
#define HPX_COMMUNICATOR_HPP

#include "../general_definitions.hpp"

#include <hpx/hpx.hpp>
#include <hpx/include/iostreams.hpp>

#include "preprocessor/mesh_metadata.hpp"

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
    HPXCommunicator(const std::string& neighborhood_data_file, const uint locality_id, const uint submesh_id);

    uint GetRankBoundaryNumber() { return this->rank_boundaries.size(); }
    HPXRankBoundary& GetRankBoundary(uint rank_boundary_id) { return this->rank_boundaries.at(rank_boundary_id); }

    void SendAll(const uint timestamp);
    hpx::future<void> ReceiveAll(const uint timestamp);

  public:
    using RankBoundaryType = HPXRankBoundary;
};

#endif
