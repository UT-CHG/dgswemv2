#ifndef SWE_DISTRIBUTED_HPP
#define SWE_DISTRIBUTED_HPP

#include <vector>

namespace SWE {
class Distributed {
  private:
    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

    uint ze_in_index;
    uint qx_in_index;
    uint qy_in_index;

    uint ze_ex_index;
    uint qx_ex_index;
    uint qy_ex_index;

  public:
    Distributed(std::vector<double>& send_buffer,
                std::vector<double>& receive_buffer,
    uint ze_in_index,
    uint qx_in_index,
    uint qy_in_index,
    uint ze_ex_index,
    uint qx_ex_index,
    uint qy_ex_index)
        : send_buffer(send_buffer),
          receive_buffer(recv_buffer),
          ze_in_index(ze_in_index),
          qx_in_index(qx_in_index),
          qy_in_index(qy_in_index),
          ze_ex_index(ze_ex_index),
          qx_ex_index(qx_ex_index),
          qy_ex_index(qy_ex_index) {}

        void GetEX(const Stepper& stepper,
                uint gp,
                const Array2D<double>& surface_normal,
                const std::vector<double>& ze_in,
                const std::vector<double>& qx_in,
                const std::vector<double>& qy_in,
                double& ze_ex,
                double& qx_ex,
                double& qy_ex){
                    ze_ex = this->receive_buffer[ze_ex_index - gp];
                    qx_ex = this->receive_buffer[qx_ex_index - gp];
                    qy_ex = this->receive_buffer[qy_ex_index - gp];
                }    

        void SetEX(
                const std::vector<double>& ze_in,
                const std::vector<double>& qx_in,
                const std::vector<double>& qy_in){
                    for(uint gp=0; gp<ze_in.size();gp++){
                        this->send_buffer[ze_in_index+gp] = ze_in[gp];
                        this->send_buffer[qx_in_index+gp] = qx_in[gp];
                        this->send_buffer[qy_in_index+gp] = qy_in[gp];
                    }
                }
};
}
#endif