#include "hpx_communicator.hpp"

HPX_REGISTER_CHANNEL(array_double);

HPXCommunicator::HPXCommunicator(const std::string& neighborhood_data_file,
                                 const uint locality_id,
                                 const uint submesh_id) {
    std::ifstream file(neighborhood_data_file);

    if (!file) {
        throw std::logic_error("Error: Unable to find distributed boundary file : " + neighborhood_data_file + '\n');
    }

    std::string line;

    uint locality_A, locality_B, submesh_A, submesh_B, n_dboubdaries;

    while (std::getline(file, line)) {
        std::stringstream neighborhood_data(line);

        neighborhood_data >> locality_A >> submesh_A >> locality_B >> submesh_B >> n_dboubdaries;

        HPXRankBoundary rank_boundary;

        std::string my_location;
        std::string neighbor_location;

        if (locality_A == locality_id && submesh_A == submesh_id) {
            my_location = std::to_string(locality_A) + "_" + std::to_string(submesh_A);
            neighbor_location = std::to_string(locality_B) + "_" + std::to_string(submesh_B);
        } else {
            my_location = std::to_string(locality_B) + "_" + std::to_string(submesh_B);
            neighbor_location = std::to_string(locality_A) + "_" + std::to_string(submesh_A);
        }

        std::string outgoing_channel_string = "channel_from_" + my_location + "_to_" + neighbor_location;
        std::string incoming_channel_string = "channel_from_" + neighbor_location + "_to_" + my_location;

        rank_boundary.outgoing = hpx::lcos::channel<array_double>(hpx::find_here());
        hpx::future<void> set_outgoing_channel = rank_boundary.outgoing.register_as(outgoing_channel_string);

        rank_boundary.incoming.connect_to(incoming_channel_string);

        rank_boundary.elements.reserve(n_dboubdaries);
        rank_boundary.bound_ids.reserve(n_dboubdaries);
        rank_boundary.p.reserve(n_dboubdaries);

        for (uint db = 0; db < n_dboubdaries; ++db) {
            DistributedBoundaryMetaData dboundary_meta_data;

            file >> dboundary_meta_data;

            file.ignore(1000, '\n');

            if (locality_A == locality_id && submesh_A == submesh_id) {
                rank_boundary.elements.push_back(dboundary_meta_data.elements.first);
                rank_boundary.bound_ids.push_back(dboundary_meta_data.bound_ids.first);
                rank_boundary.p.push_back(dboundary_meta_data.p);
            } else {
                rank_boundary.elements.push_back(dboundary_meta_data.elements.second);
                rank_boundary.bound_ids.push_back(dboundary_meta_data.bound_ids.second);
                rank_boundary.p.push_back(dboundary_meta_data.p);
            }
        }

        this->rank_boundaries.push_back(std::move(rank_boundary));

        set_outgoing_channel.get();
    }
}

void HPXCommunicator::SendPreprocAll(const uint timestamp) {
    for (auto& rank_boundary : this->rank_boundaries) {
        rank_boundary.send_preproc(timestamp);
    }
}

hpx::future<void> HPXCommunicator::ReceivePreprocAll(const uint timestamp) {
    std::vector<hpx::future<void>> receive_futures;
    receive_futures.reserve(this->rank_boundaries.size());

    for (auto& rank_boundary : this->rank_boundaries) {
        receive_futures.push_back(rank_boundary.receive_preproc(timestamp));
    }

    return hpx::when_all(receive_futures);
}

void HPXCommunicator::SendAll(const uint timestamp) {
    for (auto& rank_boundary : this->rank_boundaries) {
        rank_boundary.send(timestamp);
    }
}

hpx::future<void> HPXCommunicator::ReceiveAll(const uint timestamp) {
    std::vector<hpx::future<void>> receive_futures;
    receive_futures.reserve(this->rank_boundaries.size());

    for (auto& rank_boundary : this->rank_boundaries) {
        receive_futures.push_back(rank_boundary.receive(timestamp));
    }

    return hpx::when_all(receive_futures);
}

void HPXCommunicator::SendPostprocAll(const uint timestamp) {
    for (auto& rank_boundary : this->rank_boundaries) {
        rank_boundary.send_postproc(timestamp);
    }
}

hpx::future<void> HPXCommunicator::ReceivePostprocAll(const uint timestamp) {
    std::vector<hpx::future<void>> receive_futures;
    receive_futures.reserve(this->rank_boundaries.size());

    for (auto& rank_boundary : this->rank_boundaries) {
        receive_futures.push_back(rank_boundary.receive_postproc(timestamp));
    }

    return hpx::when_all(receive_futures);
}