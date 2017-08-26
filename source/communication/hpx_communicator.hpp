#ifndef HPX_COMMUNICATOR_HPP
#define HPX_COMMUNICATOR_HPP

#include <hpx/hpx.hpp>
#include <vector>

using vec_type = std::vector<double>;
HPX_REGISTER_CHANNEL_DECLARATION(vec_type);

struct RankInterface {
    uint locality_id;
    uint sbmsh_id;
    std::vector<uint> elements;
    std::vector<uint> face_ids;
    std::vector<uint> polynomial_order;
    std::vector<double> send_buffer;
    std::vector<double> receive_buffer;

    hpx::lcos::channel<vec_type> outgoing;
    hpx::lcos::channel<vec_type> incoming;

    inline void send(uint timestamp) { outgoing.set(send_buffer, timestamp); }

    inline hpx::future<void> receive(uint timestamp) {
        return incoming.get(timestamp)
            .then([this](hpx::future<vec_type> msg_future) { this->receive_buffer = msg_future.get(); });
    }
};

class HPXCommunicator {
  private:
    std::vector<RankInterface> distributed_interfaces;

  public:
    HPXCommunicator() = default;
    HPXCommunicator(const uint locality_id, const uint sbmsh_id, const std::string& distributed_interface_file);

    template <typename IntegrationType>
    void ResizeBuffers(IntegrationType integration, const uint num_fields);

    std::vector<double>& GetSendBufferReference(uint neighbor_id);

    std::vector<double>& GetReceiveBufferReference(uint neighbor_id);

    uint GetNumNeighbors();

    uint GetNumEdges(uint neighbor);

    std::tuple<uint, uint, uint> GetElt_FaceID_PolynomialOrder(uint neighbor, uint edge);

    void send_all(uint timestamp);

    hpx::future<void> receive_all(uint timestamp);
};

template <typename IntegrationType>
void HPXCommunicator::ResizeBuffers(IntegrationType integration, const uint num_fields) {
    for (auto& di : distributed_interfaces) {
        std::size_t sz = 0;
        for (uint i = 0; i < di.elements.size(); ++i) {
            sz += num_fields * integration.GetNumGP(di.polynomial_order[i]);
        }
        di.send_buffer.resize(sz);
        di.receive_buffer.resize(sz);
    }
}
#endif