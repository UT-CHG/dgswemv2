#ifndef DB_DATA_EXCH_HPP
#define DB_DATA_EXCH_HPP

#include "general_definitions.hpp"

class DBDataExchanger {
  public:
    const uint locality_in;
    const uint submesh_in;

    const uint locality_ex;
    const uint submesh_ex;

  private:
    std::vector<uint> offset;

    std::vector<std::vector<double>>& send_buffer;
    std::vector<std::vector<double>>& receive_buffer;

  public:
    DBDataExchanger() = default;
    DBDataExchanger(const uint locality_in,
                    const uint submesh_in,
                    const uint locality_ex,
                    const uint submesh_ex,
                    std::vector<uint>&& offset,
                    std::vector<std::vector<double>>& send_buffer,
                    std::vector<std::vector<double>>& receive_buffer);

    void SetToSendBuffer(const uint comm_type, const std::vector<double>& message);
    void GetFromReceiveBuffer(const uint comm_type, std::vector<double>& message);
};

DBDataExchanger::DBDataExchanger(const uint locality_in,
                                 const uint submesh_in,
                                 const uint locality_ex,
                                 const uint submesh_ex,
                                 std::vector<uint>&& offset,
                                 std::vector<std::vector<double>>& send_buffer,
                                 std::vector<std::vector<double>>& receive_buffer)
    : locality_in(locality_in),
      submesh_in(submesh_in),
      locality_ex(locality_ex),
      submesh_ex(submesh_ex),
      offset(std::move(offset)),
      send_buffer(send_buffer),
      receive_buffer(receive_buffer) {}

void DBDataExchanger::SetToSendBuffer(const uint comm_type, const std::vector<double>& message) {
    std::copy(message.begin(), message.end(), this->send_buffer[comm_type].begin() + this->offset[comm_type]);
}

void DBDataExchanger::GetFromReceiveBuffer(const uint comm_type, std::vector<double>& message) {
    std::copy(this->receive_buffer[comm_type].begin() + this->offset[comm_type],
              this->receive_buffer[comm_type].begin() + this->offset[comm_type] + message.size(),
              message.begin());
}

#endif
