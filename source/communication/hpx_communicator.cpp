#include "hpx_communicator.hpp"
#include "source/preprocessor/mesh_metadata.hpp"

HpxCommunicator::HpxCommunicator( const uint locality_id,
                                  const uint sbmsh_id,
                                  const std::string& distributed_interface_file) {

    std::ifstream file(distributed_interface_file.c_str());
    if ( !file ) {
        throw std::exception("Error: Unable to find distributed interface file : "
                             + distributed_interface_file + '\n');
    }

    std::string line;
    uint loc_A, loc_B, sbmsh_A, sbmsh_B, num_faces;
    while( std::getline(file, line) ) {
        std::stringstream ss(line);
        ss >> loc_A >> sbmsh_A >> loc_B >> sbmsh_B >> num_faces;
        if ( ( loc_A == locality_id && sbmsh_A == sbmsh_id ) ||
             ( loc_B == locality_id && sbmsh_B == sbmsh_id ) ) {

            uint is_first = ( loc_B == locality_id && sbmsh_B == sbmsh_id );

            RankInterface neigh_interface;

            if ( is_first ) {
                neigh_interface.locality_id = loc_B;
                neigh_interface.sbmsh_id = sbmsh_B;
            } else {
                neigh_interface.locality_id = loc_A;
                neigh_interface.sbmsh_id = sbmsh_A;
            }

            neigh_interface.elements.reserve(num_faces);
            neigh_interface.face_ids.reserve(num_faces);
            neigh_interface.polynomial_order.reserve(num_faces);

            for ( uint l = 0; l < num_faces; ++l ) {
                DistributedInterfaceMetaData dist_int;
                file >> dist_int;
                file.ignore(1000,'\n');
                if ( is_first ) {
                    neigh_interface.elements.push_back(dist_int.elements.first);
                    neigh_interface.face_id.push_back(dist_int.face_id.first);
                    neigh_interface.polynomial_order.push_back(dist_int.polynomial_order);
                } else {
                    neigh_interface.elements.push_back(dist_int.elements.second);
                    neigh_interface.face_id.push_back(dist_int.face_id.second);
                    neigh_interface.polynomial_order.push_back(dist_int.polynomial_order);
                }
            }

            distributed_interfaces.push_back( std::move(neigh_interface) );

        } else {
            for ( uint l = 0; l < num_faces; ++l ) {
                std::getline(file,line);
            }
        }
    }

    file.close();
}

const std::vector<RankInterface>& HpxCommunicator::Get_distributed_interfaces() {
    return distributed_interfaces;
}

void HpxCommunicator::resize_buffer(uint locality_id, uint sbmsh_id, std::size_t sz) {
    for ( RankInterface& rnk : distributed_interfaces ) {
        if ( rnk.locality_id == locality_id && rnk.sbmsh_id == sbmsh_id ) {
            rnk.send_buffer.resize(sz);
            rnk.receive_buffer.resize(sz);
            return;
        }
    }
}