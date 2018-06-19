#ifndef SWE_DBC_DB_DATA_EXCH_HPP
#define SWE_DBC_DB_DATA_EXCH_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/rkdg_simulation/rk_stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace DBC {
struct DBIndex {
    uint x_at_baryctr;
    uint y_at_baryctr;

    uint wet_dry;

    uint ze_in;
    uint qx_in;
    uint qy_in;

    uint ze_ex;
    uint qx_ex;
    uint qy_ex;

    uint wet_dry_postproc;

    uint ze_at_baryctr;
    uint qx_at_baryctr;
    uint qy_at_baryctr;
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

    void SetPreprocEX(const double x_at_baryctr_in, const double y_at_baryctr_in);

    void SetWetDryEX(const bool wet_in);
    void SetEX(const std::vector<double>& ze_in, const std::vector<double>& qx_in, const std::vector<double>& qy_in);

    void SetPostprocWetDryEX(const bool wet_in);
    void SetPostprocEX(const double ze_at_baryctr_in, const double qx_at_baryctr_in, const double qy_at_baryctr_in);

    void GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex);

    void GetWetDryEX(bool& wet_ex);
    void GetEX(const uint gp, double& ze_ex, double& qx_ex, double& qy_ex);

    void GetPostprocWetDryEX(bool& wet_ex);
    void GetPostprocEX(double& ze_at_baryctr_ex, double& qx_at_baryctr_ex, double& qy_at_baryctr_ex);
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

void DBDataExchanger::SetPreprocEX(const double x_at_baryctr_in, const double y_at_baryctr_in) {
    this->send_preproc_buffer[this->index.x_at_baryctr] = x_at_baryctr_in;
    this->send_preproc_buffer[this->index.y_at_baryctr] = y_at_baryctr_in;
}

void DBDataExchanger::SetWetDryEX(const bool wet_in) {
    this->send_buffer[this->index.wet_dry] = (double)wet_in;
}

void DBDataExchanger::SetEX(const std::vector<double>& ze_in,
                            const std::vector<double>& qx_in,
                            const std::vector<double>& qy_in) {
    for (uint gp = 0; gp < ze_in.size(); gp++) {
        this->send_buffer[this->index.ze_in + gp] = ze_in[gp];
        this->send_buffer[this->index.qx_in + gp] = qx_in[gp];
        this->send_buffer[this->index.qy_in + gp] = qy_in[gp];
    }
}

void DBDataExchanger::SetPostprocWetDryEX(const bool wet_in) {
    this->send_postproc_buffer[this->index.wet_dry_postproc] = (double)wet_in;
}

void DBDataExchanger::SetPostprocEX(const double ze_at_baryctr_in,
                                    const double qx_at_baryctr_in,
                                    const double qy_at_baryctr_in) {
    this->send_postproc_buffer[this->index.ze_at_baryctr] = ze_at_baryctr_in;
    this->send_postproc_buffer[this->index.qx_at_baryctr] = qx_at_baryctr_in;
    this->send_postproc_buffer[this->index.qy_at_baryctr] = qy_at_baryctr_in;
}

void DBDataExchanger::GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex) {
    x_at_baryctr_ex = this->receive_preproc_buffer[this->index.x_at_baryctr];
    y_at_baryctr_ex = this->receive_preproc_buffer[this->index.y_at_baryctr];
}

void DBDataExchanger::GetWetDryEX(bool& wet_ex) {
    wet_ex = (bool)this->receive_buffer[this->index.wet_dry];
}

void DBDataExchanger::GetEX(const uint gp, double& ze_ex, double& qx_ex, double& qy_ex) {
    ze_ex = this->receive_buffer[this->index.ze_ex - gp];
    qx_ex = this->receive_buffer[this->index.qx_ex - gp];
    qy_ex = this->receive_buffer[this->index.qy_ex - gp];
}

void DBDataExchanger::GetPostprocWetDryEX(bool& wet_ex) {
    wet_ex = (bool)this->receive_postproc_buffer[this->index.wet_dry_postproc];
}

void DBDataExchanger::GetPostprocEX(double& ze_at_baryctr_ex, double& qx_at_baryctr_ex, double& qy_at_baryctr_ex) {
    ze_at_baryctr_ex = this->receive_postproc_buffer[this->index.ze_at_baryctr];
    qx_at_baryctr_ex = this->receive_postproc_buffer[this->index.qx_at_baryctr];
    qy_at_baryctr_ex = this->receive_postproc_buffer[this->index.qy_at_baryctr];
}
}
}

#endif
