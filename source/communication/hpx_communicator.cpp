#include "hpx_communicator.hpp"
#include "preprocessor/mesh_metadata.hpp"

#include <hpx/include/iostreams.hpp>

HPX_REGISTER_CHANNEL(vec_type);

HPXCommunicator::HPXCommunicator(const uint locality_id,
                                 const uint submesh_id,
                                 const std::string& distributed_interface_file) {

    std::ifstream file(distributed_interface_file.c_str());
    if (!file) {
        throw std::logic_error("Error: Unable to find distributed interface file : " + distributed_interface_file +
                               '\n');
    }

    hpx::cout << "In communicator ctor: " << distributed_interface_file << '\n';
    std::string line;
    uint loc_A, loc_B, sbmsh_A, sbmsh_B, num_faces;
    while (std::getline(file, line)) {
        hpx::cout << "Hello " + line + "\n";
        std::stringstream ss(line);
        ss >> loc_A >> sbmsh_A >> loc_B >> sbmsh_B >> num_faces;

        uint is_first = (loc_B == locality_id && sbmsh_B == submesh_id);

        RankInterface neigh_interface;
        neigh_interface.outgoing = hpx::lcos::channel<vec_type>(hpx::find_here());

        std::string in_location;
        std::string out_location;

        if (is_first) {
            neigh_interface.locality_id = loc_B;
            neigh_interface.submesh_id = sbmsh_B;

            in_location = std::to_string(loc_A) + "_" + std::to_string(sbmsh_A);
            out_location = std::to_string(loc_B) + "_" + std::to_string(sbmsh_B);
        } else {
            neigh_interface.locality_id = loc_A;
            neigh_interface.submesh_id = sbmsh_A;

            in_location = std::to_string(loc_B) + "_" + std::to_string(sbmsh_B);
            out_location = std::to_string(loc_A) + "_" + std::to_string(sbmsh_A);
        }

        std::string outgoing_channel_string = "channel_from_" + in_location + "_to_" + out_location;
        std::string incoming_channel_string = "channel_from_" + out_location + "_to_" + in_location;

        hpx::cout << "Setting up outgoing channel: " << outgoing_channel_string << hpx::endl;
        hpx::cout << "Setting up incoming channel: " << incoming_channel_string << hpx::endl;

        hpx::future<void> set_out = neigh_interface.outgoing.register_as(outgoing_channel_string);
        neigh_interface.incoming.connect_to(incoming_channel_string);

        neigh_interface.elements.reserve(num_faces);
        neigh_interface.face_ids.reserve(num_faces);
        neigh_interface.polynomial_order.reserve(num_faces);

        for (uint l = 0; l < num_faces; ++l) {
            DistributedInterfaceMetaData dist_int;
            file >> dist_int;
            file.ignore(1000, '\n');
            if (is_first) {
                neigh_interface.elements.push_back(dist_int.elements.first);
                neigh_interface.face_ids.push_back(dist_int.face_id.first);
                neigh_interface.polynomial_order.push_back(dist_int.polynomial_order);
            } else {
                neigh_interface.elements.push_back(dist_int.elements.second);
                neigh_interface.face_ids.push_back(dist_int.face_id.second);
                neigh_interface.polynomial_order.push_back(dist_int.polynomial_order);
            }
        }

        distributed_interfaces.push_back(std::move(neigh_interface));

        set_out.get();
    }

    file.close();
}

std::vector<double>& HPXCommunicator::GetSendBufferReference(uint neighbor_id) {
    return distributed_interfaces.at(neighbor_id).send_buffer;
}

std::vector<double>& HPXCommunicator::GetReceiveBufferReference(uint neighbor_id) {
    return distributed_interfaces.at(neighbor_id).receive_buffer;
}

uint HPXCommunicator::GetNumNeighbors() { return distributed_interfaces.size(); }

uint HPXCommunicator::GetNumEdges(uint neighbor) { return distributed_interfaces.at(neighbor).elements.size(); }

std::tuple<uint, uint, uint> HPXCommunicator::GetElt_FaceID_PolynomialOrder(uint neighbor, uint edge) {
    return std::make_tuple(distributed_interfaces[neighbor].elements[edge],
                           distributed_interfaces[neighbor].face_ids[edge],
                           distributed_interfaces[neighbor].polynomial_order[edge]);
}

void HPXCommunicator::send_all(uint timestamp) {
    for (auto&& rnk : distributed_interfaces) {
        rnk.send(timestamp);
    }
}

hpx::future<void> HPXCommunicator::receive_all(uint timestamp) {
    std::vector<hpx::future<void>> receive_futures;
    receive_futures.reserve(distributed_interfaces.size());

    for (auto&& rnk : distributed_interfaces) {
        receive_futures.push_back(rnk.receive(timestamp));
    }

    return hpx::when_all(receive_futures);
}