#include "hpx_communicator.hpp"

HPX_REGISTER_CHANNEL(array_double);

HPXCommunicator::HPXCommunicator(const DistributedBoundaryMetaData& db_data) {
    for (auto& rb_meta_data : db_data.rank_boundary_data) {
        HPXRankBoundary rank_boundary;

        rank_boundary.db_data = rb_meta_data;

        std::string my_location;
        std::string neighbor_location;

        my_location       = std::to_string(rb_meta_data.locality_in) + "_" + std::to_string(rb_meta_data.submesh_in);
        neighbor_location = std::to_string(rb_meta_data.locality_ex) + "_" + std::to_string(rb_meta_data.submesh_ex);

        std::string outgoing_channel_string = "channel_from_" + my_location + "_to_" + neighbor_location;
        std::string incoming_channel_string = "channel_from_" + neighbor_location + "_to_" + my_location;

        rank_boundary.outgoing                 = hpx::lcos::channel<array_double>(hpx::find_here());
        hpx::future<void> set_outgoing_channel = rank_boundary.outgoing.register_as(outgoing_channel_string);

        rank_boundary.incoming.connect_to(incoming_channel_string);

        this->rank_boundaries.push_back(std::move(rank_boundary));

        set_outgoing_channel.get();
    }
}

void HPXCommunicator::SendAll(const uint comm_type, const uint timestamp) {
    for (auto& rank_boundary : this->rank_boundaries) {
        rank_boundary.outgoing.set(rank_boundary.send_buffer[comm_type], timestamp);
    }
}

hpx::future<void> HPXCommunicator::ReceiveAll(const uint comm_type, const uint timestamp) {
    std::vector<hpx::future<void>> receive_futures;
    receive_futures.reserve(this->rank_boundaries.size());

    for (auto& rank_boundary : this->rank_boundaries) {
        receive_futures.push_back(rank_boundary.incoming.get(timestamp).then(
            [&rank_boundary, comm_type](hpx::future<array_double> msg_future) {
                rank_boundary.receive_buffer[comm_type] = msg_future.get();
            }));
    }

    return hpx::when_all(receive_futures);
}