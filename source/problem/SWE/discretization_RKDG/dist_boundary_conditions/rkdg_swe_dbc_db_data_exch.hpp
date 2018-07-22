#ifndef RKDG_SWE_DBC_DB_DATA_EXCH_HPP
#define RKDG_SWE_DBC_DB_DATA_EXCH_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace DBC {
struct DBIndex {
    uint x_at_baryctr;
    uint y_at_baryctr;

    uint wet_dry;
    uint q_in;
    uint q_ex;

    uint wet_dry_postproc;
    uint q_at_baryctr;
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
    void SetEX(const DynMatrix<double>& q_in);

    void SetPostprocWetDryEX(const bool wet_in);
    void SetPostprocEX(const StatVector<double, SWE::n_variables>& q_at_baryctr_in);

    void GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex);

    void GetWetDryEX(bool& wet_ex);
    void GetEX(const uint gp, DynVector<double>& q_ex);

    void GetPostprocWetDryEX(bool& wet_ex);
    void GetPostprocEX(StatVector<double, SWE::n_variables>& q_at_baryctr_ex);
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

void DBDataExchanger::SetEX(const DynMatrix<double>& q_in) {
    for (uint gp = 0; gp < columns(q_in); gp++) {
        for (uint var = 0; var < SWE::n_variables; var++) {
            this->send_buffer[this->index.q_in + SWE::n_variables * gp + var] = q_in(var, gp);
        }
    }
}

void DBDataExchanger::SetPostprocWetDryEX(const bool wet_in) {
    this->send_postproc_buffer[this->index.wet_dry_postproc] = (double)wet_in;
}

void DBDataExchanger::SetPostprocEX(const StatVector<double, SWE::n_variables>& q_at_baryctr_in) {
    for (uint var = 0; var < SWE::n_variables; var++) {
        this->send_postproc_buffer[this->index.q_at_baryctr + var] = q_at_baryctr_in[var];
    }
}

void DBDataExchanger::GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex) {
    x_at_baryctr_ex = this->receive_preproc_buffer[this->index.x_at_baryctr];
    y_at_baryctr_ex = this->receive_preproc_buffer[this->index.y_at_baryctr];
}

void DBDataExchanger::GetWetDryEX(bool& wet_ex) {
    wet_ex = (bool)this->receive_buffer[this->index.wet_dry];
}

void DBDataExchanger::GetEX(const uint gp, DynVector<double>& q_ex) {
    for (uint var = 0; var < SWE::n_variables; var++) {
        q_ex[SWE::n_variables - var - 1] = this->receive_buffer[this->index.q_ex - SWE::n_variables * gp - var];
    }
}

void DBDataExchanger::GetPostprocWetDryEX(bool& wet_ex) {
    wet_ex = (bool)this->receive_postproc_buffer[this->index.wet_dry_postproc];
}

void DBDataExchanger::GetPostprocEX(StatVector<double, SWE::n_variables>& q_at_baryctr_ex) {
    for (uint var = 0; var < SWE::n_variables; var++) {
        q_at_baryctr_ex[var] = this->receive_postproc_buffer[this->index.q_at_baryctr + var];
    }
}
}
}
}

#endif
