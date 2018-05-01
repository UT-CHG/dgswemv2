#ifndef SWE_DBC_DISTRIBUTED_HPP
#define SWE_DBC_DISTRIBUTED_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace DBC {
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
                const uint bath_at_baryctr_index);

    template <typename DistributedBoundaryType>
    void ComputeFlux(const Stepper& stepper, DistributedBoundaryType& dbound);

    void SetPreprocEX(const double x_at_baryctr_in, const double y_at_baryctr_in);
    void SetWetDryEX(const bool wet_in);
    void SetEX(const std::vector<double>& ze_in, const std::vector<double>& qx_in, const std::vector<double>& qy_in);
    void SetPostprocEX(const double ze_at_baryctr_in,
                       const double qx_at_baryctr_in,
                       const double qy_at_baryctr_in,
                       const double bath_at_baryctr_in);

    void GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex);
    void GetWetDryEX(bool& wet_ex);
    void GetEX(const Stepper& stepper,
               const uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex);
    void GetPostprocEX(double& ze_at_baryctr_ex,
                       double& qx_at_baryctr_ex,
                       double& qy_at_baryctr_ex,
                       double& bath_at_baryctr_ex);
};

Distributed::Distributed(std::vector<double>& send_preproc_buffer,
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

template <typename DistributedBoundaryType>
void Distributed::ComputeFlux(const Stepper& stepper, DistributedBoundaryType& dbound) {
    bool wet_in = dbound.data.wet_dry_state.wet;
    bool wet_ex;
    this->GetWetDryEX(wet_ex);

    const uint stage = stepper.GetStage();

    auto& state    = dbound.data.state[stage];
    auto& boundary = dbound.data.boundary[dbound.bound_id];
    auto& sp_at_gp = dbound.data.spherical_projection.sp_at_gp_boundary[dbound.bound_id];

    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
        this->GetEX(stepper,
                    gp,
                    dbound.surface_normal,
                    boundary.ze_at_gp,
                    boundary.qx_at_gp,
                    boundary.qy_at_gp,
                    ze_ex,
                    qx_ex,
                    qy_ex);

        LLF_flux(Global::g,
                 boundary.ze_at_gp[gp],
                 ze_ex,
                 boundary.qx_at_gp[gp],
                 qx_ex,
                 boundary.qy_at_gp[gp],
                 qy_ex,
                 boundary.bath_at_gp[gp],
                 sp_at_gp[gp],
                 dbound.surface_normal[gp],
                 boundary.ze_numerical_flux_at_gp[gp],
                 boundary.qx_numerical_flux_at_gp[gp],
                 boundary.qy_numerical_flux_at_gp[gp]);
    }

    // compute net volume flux out of IN/EX elements
    double net_volume_flux_in = 0;

    net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);

    if (net_volume_flux_in > 0) {
        if (!wet_in) {  // water flowing from dry IN element
            // Zero flux on IN element side
            std::fill(boundary.ze_numerical_flux_at_gp.begin(), boundary.ze_numerical_flux_at_gp.end(), 0.0);
            std::fill(boundary.qx_numerical_flux_at_gp.begin(), boundary.qx_numerical_flux_at_gp.end(), 0.0);
            std::fill(boundary.qy_numerical_flux_at_gp.begin(), boundary.qy_numerical_flux_at_gp.end(), 0.0);

            net_volume_flux_in = 0;
        } else if (!wet_ex) {  // water flowing to dry EX element
            net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);
        }
    } else if (net_volume_flux_in < 0) {
        if (!wet_ex) {  // water flowing from dry EX element
            // Reflective Boundary on IN element side
            SWE::BC::Land land_boundary;

            land_boundary.ComputeFlux(stepper,
                                      dbound.surface_normal,
                                      sp_at_gp,
                                      boundary.bath_at_gp,
                                      boundary.ze_at_gp,
                                      boundary.qx_at_gp,
                                      boundary.qy_at_gp,
                                      boundary.ze_numerical_flux_at_gp,
                                      boundary.qx_numerical_flux_at_gp,
                                      boundary.qy_numerical_flux_at_gp);

            net_volume_flux_in = 0;
        } else if (!wet_in) {  // water flowing to dry IN element
            for (uint gp = 0; gp < dbound.data.get_ngp_boundary(dbound.bound_id); ++gp) {
                this->GetEX(stepper,
                            gp,
                            dbound.surface_normal,
                            boundary.ze_at_gp,
                            boundary.qx_at_gp,
                            boundary.qy_at_gp,
                            ze_ex,
                            qx_ex,
                            qy_ex);

                LLF_flux(0.0,
                         boundary.ze_at_gp[gp],
                         ze_ex,
                         boundary.qx_at_gp[gp],
                         qx_ex,
                         boundary.qy_at_gp[gp],
                         qy_ex,
                         boundary.bath_at_gp[gp],
                         sp_at_gp[gp],
                         dbound.surface_normal[gp],
                         boundary.ze_numerical_flux_at_gp[gp],
                         boundary.qx_numerical_flux_at_gp[gp],
                         boundary.qy_numerical_flux_at_gp[gp]);
            }

            net_volume_flux_in = dbound.Integration(boundary.ze_numerical_flux_at_gp);
        }
    }
}

void Distributed::SetPreprocEX(const double x_at_baryctr_in, const double y_at_baryctr_in) {
    this->send_preproc_buffer[x_at_baryctr_index] = x_at_baryctr_in;
    this->send_preproc_buffer[y_at_baryctr_index] = y_at_baryctr_in;
}

void Distributed::SetWetDryEX(const bool wet_in) {
    this->send_buffer[wet_dry_index] = (double)wet_in;
}

void Distributed::SetEX(const std::vector<double>& ze_in,
                        const std::vector<double>& qx_in,
                        const std::vector<double>& qy_in) {
    for (uint gp = 0; gp < ze_in.size(); gp++) {
        this->send_buffer[ze_in_index + gp] = ze_in[gp];
        this->send_buffer[qx_in_index + gp] = qx_in[gp];
        this->send_buffer[qy_in_index + gp] = qy_in[gp];
    }
}

void Distributed::SetPostprocEX(const double ze_at_baryctr_in,
                                const double qx_at_baryctr_in,
                                const double qy_at_baryctr_in,
                                const double bath_at_baryctr_in) {
    this->send_postproc_buffer[ze_at_baryctr_index]   = ze_at_baryctr_in;
    this->send_postproc_buffer[qx_at_baryctr_index]   = qx_at_baryctr_in;
    this->send_postproc_buffer[qy_at_baryctr_index]   = qy_at_baryctr_in;
    this->send_postproc_buffer[bath_at_baryctr_index] = bath_at_baryctr_in;
}

void Distributed::GetPreprocEX(double& x_at_baryctr_ex, double& y_at_baryctr_ex) {
    x_at_baryctr_ex = this->receive_preproc_buffer[x_at_baryctr_index];
    y_at_baryctr_ex = this->receive_preproc_buffer[y_at_baryctr_index];
}

void Distributed::GetWetDryEX(bool& wet_ex) {
    wet_ex = (bool)this->receive_buffer[wet_dry_index];
}

void Distributed::GetEX(const Stepper& stepper,
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

void Distributed::GetPostprocEX(double& ze_at_baryctr_ex,
                                double& qx_at_baryctr_ex,
                                double& qy_at_baryctr_ex,
                                double& bath_at_baryctr_ex) {
    ze_at_baryctr_ex   = this->receive_postproc_buffer[ze_at_baryctr_index];
    qx_at_baryctr_ex   = this->receive_postproc_buffer[qx_at_baryctr_index];
    qy_at_baryctr_ex   = this->receive_postproc_buffer[qy_at_baryctr_index];
    bath_at_baryctr_ex = this->receive_postproc_buffer[bath_at_baryctr_index];
}
}
}

#endif