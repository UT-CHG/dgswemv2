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
    // std::unordered_map<uint> distributed_interface2rank;

  public:
    HPXCommunicator() = default;
    HPXCommunicator(const uint locality_id, const uint sbmsh_id, const std::string& distributed_interface_file);

    void resize_buffer(uint locality_id, uint sbmsh_id, std::size_t sz);

    inline const std::vector<RankInterface>& Get_distributed_interfaces() const { return distributed_interfaces; }

    void send_all(uint timestamp);

    hpx::future<void> receive_all(uint timestamp);
};
#endif