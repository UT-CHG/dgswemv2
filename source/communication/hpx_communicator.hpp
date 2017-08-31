#ifndef HPX_COMMUNICATOR_HPP
#define HPX_COMMUNICATOR_HPP

#include <hpx/hpx.hpp>
#include <vector>

using array_double = std::vector<double>;
HPX_REGISTER_CHANNEL_DECLARATION(array_double);

struct HPXRankBoundary {
    uint neighbor_locality_id;
    uint neighbor_submesh_id;

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
    HPXCommunicator(const std::string&, const uint, const uint);

    void SendAll(uint);
    hpx::future<void> ReceiveAll(uint);

    uint GetRankBoundaryNumber() { return rank_boundaries.size(); }

    HPXRankBoundary& GetRankBoundary(uint rank_boundary_id) { return rank_boundaries.at(rank_boundary_id); }
    /*
        template <typename IntegrationType>
        void ResizeBuffers(IntegrationType integration, const uint num_fields);

        std::vector<double>& GetSendBufferReference(uint neighbor_id);

        std::vector<double>& GetReceiveBufferReference(uint neighbor_id);

        uint GetNumNeighbors();

        uint GetNumEdges(uint neighbor);

        std::tuple<uint, uint, uint> GetElt_FaceID_PolynomialOrder(uint neighbor, uint edge);
    */
    using RankBoundaryType = HPXRankBoundary;
};
/*
template <typename IntegrationType>
void HPXCommunicator::ResizeBuffers(IntegrationType integration, const uint num_fields) {
    for (auto& di : rank_boundaries) {
        std::size_t sz = 0;
        for (uint i = 0; i < di.elements.size(); ++i) {
            sz += num_fields * integration.GetNumGP(di.p[i]);
        }
        di.send_buffer.resize(sz);
        di.receive_buffer.resize(sz);
    }
}
*/
#endif