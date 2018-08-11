#ifndef RKDG_SWE_DBC_DB_DATA_EXCH_HPP
#define RKDG_SWE_DBC_DB_DATA_EXCH_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace DBC {
class DBDataExchanger {
  private:
    std::vector<uint> offset;

    std::vector<std::vector<double>>& send_buffer;
    std::vector<std::vector<double>>& receive_buffer;

  public:
    DBDataExchanger() = default;
    DBDataExchanger(const std::vector<uint>& offset,
                    std::vector<std::vector<double>>& send_buffer,
                    std::vector<std::vector<double>>& receive_buffer);

    void SetToSendBuffer(const uint comm_type, const std::vector<double>& message);
    void GetFromReceiveBuffer(const uint comm_type, std::vector<double>& message);
};

DBDataExchanger::DBDataExchanger(const std::vector<uint>& offset,
                                 std::vector<std::vector<double>>& send_buffer,
                                 std::vector<std::vector<double>>& receive_buffer)
    : offset(offset), send_buffer(send_buffer), receive_buffer(receive_buffer) {}

void DBDataExchanger::SetToSendBuffer(const uint comm_type, const std::vector<double>& message) {
    std::copy(message.begin(), message.end(), this->send_buffer[comm_type].begin() + this->offset[comm_type]);
}

void DBDataExchanger::GetFromReceiveBuffer(const uint comm_type, std::vector<double>& message) {
    std::copy(this->receive_buffer[comm_type].begin() + this->offset[comm_type],
              this->receive_buffer[comm_type].begin() + this->offset[comm_type] + message.size(),
              message.begin());
}
}
}
}

#endif
