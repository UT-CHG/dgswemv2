#ifndef HPX_COMMUNICATOR_HPP
#define HPX_COMMUNICATOR_HPP

struct RankInterface {
    uint locality_id;
    uint sbmsh_id;
    std::vector<uint> elements;
    std::vector<uint> face_ids;
    std::vector<uint> polynomial_order;
    std::vector<double> send_buffer;
    std::vector<double> receive_buffer;
};

class HpxCommunicator {
private:
    std::vector<RankInterface> distributed_interfaces;
    //std::unordered_map<uint> distributed_interface2rank;

public:
    HpxCommunicator( const uint locality_id,
                     const uint sbmsh_id,
                     const std::string& distributed_interface_file );

    void resize_buffer(uint locality_id, uint sbmsh_id, std::size_t sz);

    inline const std::vector<RankInterface>& Get_distributed_interfaces() {
        return distributed_interfaces;
    }

    inline std::vector<double>::iterator Get_send_buffer_iterator(uint neigh, uint indx) {
        assert( indx < distributed_interfaces.at(neigh).send_buffer.size() );
        return distributed_interfaces.at(neigh).send_buffer.begin() + indx;
    }

    inline std::vector<double>::reverse_iterator Get_receive_buffer_iterator(uint neigh, uint indx) {
        assert( indx < distributed_interfaces.at(neigh).receive_buffer.size() );
        return distributed_interfaces.at(neigh).receive_buffer.rbegin() + indx;
    }

    hpx::future<void> send_all();

    hpx::future<void> receive_all();

};
#endif