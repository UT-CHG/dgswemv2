#ifndef SWE_BC_DISTRIBUTED_HPP
#define SWE_BC_DISTRIBUTED_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"

namespace SWE {
class Distributed {
  private:
    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

    std::vector<double>& send_postproc_buffer;
    std::vector<double>& receive_postproc_buffer;

    uint ze_in_index;
    uint qx_in_index;
    uint qy_in_index;

    uint ze_ex_index;
    uint qx_ex_index;
    uint qy_ex_index;

    uint ze_at_baryctr_in_index;
    uint qx_at_baryctr_in_index;
    uint qy_at_baryctr_in_index;
    uint bath_at_baryctr_in_index;

    uint ze_at_baryctr_ex_index;
    uint qx_at_baryctr_ex_index;
    uint qy_at_baryctr_ex_index;
    uint bath_at_baryctr_ex_index;

  public:
    Distributed(std::vector<double>& send_buffer,
                std::vector<double>& receive_buffer,
                std::vector<double>& send_postproc_buffer,
                std::vector<double>& receive_postproc_buffer,                
                const uint ze_in_index,
                const uint qx_in_index,
                const uint qy_in_index,
                const uint ze_ex_index,
                const uint qx_ex_index,
                const uint qy_ex_index,
                const uint ze_at_baryctr_in_index,
                const uint qx_at_baryctr_in_index,
                const uint qy_at_baryctr_in_index,
                const uint bath_at_baryctr_in_index,
                const uint ze_at_baryctr_ex_index,
                const uint qx_at_baryctr_ex_index,
                const uint qy_at_baryctr_ex_index,
                const uint bath_at_baryctr_ex_index)
        : send_buffer(send_buffer),
          receive_buffer(receive_buffer),
          send_postproc_buffer(send_postproc_buffer),
          receive_postproc_buffer(receive_postproc_buffer),
          ze_in_index(ze_in_index),
          qx_in_index(qx_in_index),
          qy_in_index(qy_in_index),
          ze_ex_index(ze_ex_index),
          qx_ex_index(qx_ex_index),
          qy_ex_index(qy_ex_index), 
          ze_at_baryctr_in_index(ze_at_baryctr_in_index),
          qx_at_baryctr_in_index(qx_at_baryctr_in_index),
          qy_at_baryctr_in_index(qy_at_baryctr_in_index),
          bath_at_baryctr_in_index(bath_at_baryctr_in_index),
          ze_at_baryctr_ex_index(ze_at_baryctr_ex_index),
          qx_at_baryctr_ex_index(qx_at_baryctr_ex_index),
          qy_at_baryctr_ex_index(qy_at_baryctr_ex_index),
          bath_at_baryctr_ex_index(bath_at_baryctr_ex_index) {}

    void SetEX(const std::vector<double>& ze_in, const std::vector<double>& qx_in, const std::vector<double>& qy_in) {
        for (uint gp = 0; gp < ze_in.size(); gp++) {
            this->send_buffer[ze_in_index + gp] = ze_in[gp];
            this->send_buffer[qx_in_index + gp] = qx_in[gp];
            this->send_buffer[qy_in_index + gp] = qy_in[gp];
        }
    }

    void GetEX(const Stepper& stepper,
               const uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex) {
        ze_ex = this->receive_buffer[ze_ex_index - gp];
        qx_ex = this->receive_buffer[qx_ex_index - gp];
        qy_ex = this->receive_buffer[qy_ex_index - gp];
    }

    void SetPostprocEX(const double ze_at_baryctr_in, const double qx_at_baryctr_in, const double qy_at_baryctr_in, const double bath_at_baryctr_in) {
            this->send_postproc_buffer[ze_at_baryctr_in_index] = ze_at_baryctr_in;
            this->send_postproc_buffer[qx_at_baryctr_in_index] = qx_at_baryctr_in;
            this->send_postproc_buffer[qy_at_baryctr_in_index] = qy_at_baryctr_in;
            this->send_postproc_buffer[bath_at_baryctr_in_index] = bath_at_baryctr_in;
    }

    void GetPostprocEX(double& ze_at_baryctr_ex, double& qx_at_baryctr_ex, double& qy_at_baryctr_ex, double& bath_at_baryctr_ex) {
        ze_at_baryctr_ex = this->receive_postproc_buffer[ze_at_baryctr_ex_index];
        qx_at_baryctr_ex = this->receive_postproc_buffer[qx_at_baryctr_ex_index];
        qy_at_baryctr_ex = this->receive_postproc_buffer[qy_at_baryctr_ex_index];
        bath_at_baryctr_ex = this->receive_postproc_buffer[bath_at_baryctr_ex_index];
    }
};
}

#endif