#ifndef SWE_BC_DISTRIBUTED_HPP
#define SWE_BC_DISTRIBUTED_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"

namespace SWE {
class Distributed {
  private:
    std::vector<double>& send_preproc_buffer;
    std::vector<double>& receive_preproc_buffer;

    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

    std::vector<double>& send_postproc_buffer;
    std::vector<double>& receive_postproc_buffer;

    uint x_at_baryctr_index;
    uint y_at_baryctr_index;

    uint wet_dry_index;

    uint ze_in_index;
    uint qx_in_index;
    uint qy_in_index;

    uint ze_ex_index;
    uint qx_ex_index;
    uint qy_ex_index;

    uint ze_at_baryctr_index;
    uint qx_at_baryctr_index;
    uint qy_at_baryctr_index;
    uint bath_at_baryctr_index;

  public:
    Distributed(std::vector<double>& send_preproc_buffer,
                std::vector<double>& receive_preproc_buffer,
                std::vector<double>& send_buffer,
                std::vector<double>& receive_buffer,
                std::vector<double>& send_postproc_buffer,
                std::vector<double>& receive_postproc_buffer,
                const uint x_at_baryctr_index,
                const uint y_at_baryctr_index,
                const uint wet_dry_index,
                const uint ze_in_index,
                const uint qx_in_index,
                const uint qy_in_index,
                const uint ze_ex_index,
                const uint qx_ex_index,
                const uint qy_ex_index,
                const uint ze_at_baryctr_index,
                const uint qx_at_baryctr_index,
                const uint qy_at_baryctr_index,
                const uint bath_at_baryctr_index)
        : send_preproc_buffer(send_preproc_buffer),
          receive_preproc_buffer(receive_preproc_buffer),
          send_buffer(send_buffer),
          receive_buffer(receive_buffer),
          send_postproc_buffer(send_postproc_buffer),
          receive_postproc_buffer(receive_postproc_buffer),
          x_at_baryctr_index(x_at_baryctr_index),
          y_at_baryctr_index(y_at_baryctr_index),
          wet_dry_index(wet_dry_index),
          ze_in_index(ze_in_index),
          qx_in_index(qx_in_index),
          qy_in_index(qy_in_index),
          ze_ex_index(ze_ex_index),
          qx_ex_index(qx_ex_index),
          qy_ex_index(qy_ex_index),
          ze_at_baryctr_index(ze_at_baryctr_index),
          qx_at_baryctr_index(qx_at_baryctr_index),
          qy_at_baryctr_index(qy_at_baryctr_index),
          bath_at_baryctr_index(bath_at_baryctr_index) {}

    void SetPreprocEX(const double x_at_baryctr_in, const double y_at_baryctr_in) {
        this->send_preproc_buffer[x_at_baryctr_index] = x_at_baryctr_in;
        this->send_preproc_buffer[y_at_baryctr_index] = y_at_baryctr_in;
    }

    void GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex) {
        x_at_baryctr_ex = this->receive_preproc_buffer[x_at_baryctr_index];
        y_at_baryctr_ex = this->receive_preproc_buffer[y_at_baryctr_index];
    }

    void SetWetDryEX(const bool wet_in) { this->send_buffer[wet_dry_index] = (double)wet_in; }

    void GetWetDryEX(bool& wet_ex) { wet_ex = (bool)this->receive_buffer[wet_dry_index]; }

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

    void SetPostprocEX(const double ze_at_baryctr_in,
                       const double qx_at_baryctr_in,
                       const double qy_at_baryctr_in,
                       const double bath_at_baryctr_in) {
        this->send_postproc_buffer[ze_at_baryctr_index] = ze_at_baryctr_in;
        this->send_postproc_buffer[qx_at_baryctr_index] = qx_at_baryctr_in;
        this->send_postproc_buffer[qy_at_baryctr_index] = qy_at_baryctr_in;
        this->send_postproc_buffer[bath_at_baryctr_index] = bath_at_baryctr_in;
    }

    void GetPostprocEX(double& ze_at_baryctr_ex,
                       double& qx_at_baryctr_ex,
                       double& qy_at_baryctr_ex,
                       double& bath_at_baryctr_ex) {
        ze_at_baryctr_ex = this->receive_postproc_buffer[ze_at_baryctr_index];
        qx_at_baryctr_ex = this->receive_postproc_buffer[qx_at_baryctr_index];
        qy_at_baryctr_ex = this->receive_postproc_buffer[qy_at_baryctr_index];
        bath_at_baryctr_ex = this->receive_postproc_buffer[bath_at_baryctr_index];
    }
};
}

#endif