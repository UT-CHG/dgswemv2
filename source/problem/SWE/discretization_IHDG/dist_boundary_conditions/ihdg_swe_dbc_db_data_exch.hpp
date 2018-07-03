#ifndef IHDG_SWE_DBC_DB_DATA_EXCH_HPP
#define IHDG_SWE_DBC_DB_DATA_EXCH_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"

namespace SWE {
namespace IHDG {
namespace DBC {
struct DBIndex {
    uint ze_in;
    uint qx_in;
    uint qy_in;
    uint ze_flux_dot_n_in;
    uint qx_flux_dot_n_in;
    uint qy_flux_dot_n_in;

    uint ze_ex;
    uint qx_ex;
    uint qy_ex;
    uint ze_flux_dot_n_ex;
    uint qx_flux_dot_n_ex;
    uint qy_flux_dot_n_ex;
};

class DBDataExchanger {
  private:
    DBIndex index;

    std::vector<double>& send_preproc_buffer;
    std::vector<double>& receive_preproc_buffer;

    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

    std::vector<double>& send_postproc_buffer;
    std::vector<double>& receive_postproc_buffer;

  public:
    DBDataExchanger() = default;
    DBDataExchanger(const DBIndex& index,
                    std::vector<double>& send_preproc_buffer,
                    std::vector<double>& receive_preproc_buffer,
                    std::vector<double>& send_buffer,
                    std::vector<double>& receive_buffer,
                    std::vector<double>& send_postproc_buffer,
                    std::vector<double>& receive_postproc_buffer);

    void SetEX(const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               const std::vector<double>& ze_flux_dot_n_in,
               const std::vector<double>& qx_flux_dot_n_in,
               const std::vector<double>& qy_flux_dot_n_in);

    void GetEX(const uint gp,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex,
               double& ze_flux_dot_n_ex,
               double& qx_flux_dot_n_ex,
               double& qy_flux_dot_n_ex);
};

DBDataExchanger::DBDataExchanger(const DBIndex& index,
                                 std::vector<double>& send_preproc_buffer,
                                 std::vector<double>& receive_preproc_buffer,
                                 std::vector<double>& send_buffer,
                                 std::vector<double>& receive_buffer,
                                 std::vector<double>& send_postproc_buffer,
                                 std::vector<double>& receive_postproc_buffer)
    : index(index),
      send_preproc_buffer(send_preproc_buffer),
      receive_preproc_buffer(receive_preproc_buffer),
      send_buffer(send_buffer),
      receive_buffer(receive_buffer),
      send_postproc_buffer(send_postproc_buffer),
      receive_postproc_buffer(receive_postproc_buffer) {}

void DBDataExchanger::SetEX(const std::vector<double>& ze_in,
                            const std::vector<double>& qx_in,
                            const std::vector<double>& qy_in,
                            const std::vector<double>& ze_flux_dot_n_in,
                            const std::vector<double>& qx_flux_dot_n_in,
                            const std::vector<double>& qy_flux_dot_n_in) {
    for (uint gp = 0; gp < ze_in.size(); gp++) {
        this->send_buffer[this->index.ze_in + gp]            = ze_in[gp];
        this->send_buffer[this->index.qx_in + gp]            = qx_in[gp];
        this->send_buffer[this->index.qy_in + gp]            = qy_in[gp];
        this->send_buffer[this->index.ze_flux_dot_n_in + gp] = ze_flux_dot_n_in[gp];
        this->send_buffer[this->index.qx_flux_dot_n_in + gp] = qx_flux_dot_n_in[gp];
        this->send_buffer[this->index.qy_flux_dot_n_in + gp] = qy_flux_dot_n_in[gp];
    }
}

void DBDataExchanger::GetEX(const uint gp,
                            double& ze_ex,
                            double& qx_ex,
                            double& qy_ex,
                            double& ze_flux_dot_n_ex,
                            double& qx_flux_dot_n_ex,
                            double& qy_flux_dot_n_ex) {
    ze_ex            = this->receive_buffer[this->index.ze_ex - gp];
    qx_ex            = this->receive_buffer[this->index.qx_ex - gp];
    qy_ex            = this->receive_buffer[this->index.qy_ex - gp];
    ze_flux_dot_n_ex = this->receive_buffer[this->index.ze_flux_dot_n_ex - gp];
    qx_flux_dot_n_ex = this->receive_buffer[this->index.qx_flux_dot_n_ex - gp];
    qy_flux_dot_n_ex = this->receive_buffer[this->index.qy_flux_dot_n_ex - gp];
}
}
}
}

#endif
