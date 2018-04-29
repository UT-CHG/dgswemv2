#ifndef SWE_BC_LAND_HPP
#define SWE_BC_LAND_HPP

#include "../../../general_definitions.hpp"
#include "../../../simulation/stepper.hpp"
#include "../numerical_fluxes/swe_numerical_fluxes.hpp"

namespace SWE {
namespace BC {
class Land {
  public:
    void ComputeFlux(const Stepper& stepper,
                     const Array2D<double>& surface_normal,
                     const std::vector<double>& sp_in,
                     const std::vector<double>& bath_in,
                     const std::vector<double>& ze_in,
                     const std::vector<double>& qx_in,
                     const std::vector<double>& qy_in,
                     std::vector<double>& ze_numerical_flux,
                     std::vector<double>& qx_numerical_flux,
                     std::vector<double>& qy_numerical_flux);

    void GetEX(const Stepper& stepper,
               const uint gp,
               const Array2D<double>& surface_normal,
               const std::vector<double>& ze_in,
               const std::vector<double>& qx_in,
               const std::vector<double>& qy_in,
               double& ze_ex,
               double& qx_ex,
               double& qy_ex);
};

void Land::ComputeFlux(const Stepper& stepper,
                       const Array2D<double>& surface_normal,
                       const std::vector<double>& sp_in,
                       const std::vector<double>& bath_in,
                       const std::vector<double>& ze_in,
                       const std::vector<double>& qx_in,
                       const std::vector<double>& qy_in,
                       std::vector<double>& ze_numerical_flux,
                       std::vector<double>& qx_numerical_flux,
                       std::vector<double>& qy_numerical_flux) {
    double ze_ex, qx_ex, qy_ex;
    for (uint gp = 0; gp < ze_in.size(); ++gp) {
        this->GetEX(stepper, gp, surface_normal, ze_in, qx_in, qy_in, ze_ex, qx_ex, qy_ex);

        LLF_flux(Global::g,
                 ze_in[gp],
                 ze_ex,
                 qx_in[gp],
                 qx_ex,
                 qy_in[gp],
                 qy_ex,
                 bath_in[gp],
                 sp_in[gp],
                 surface_normal[gp],
                 ze_numerical_flux[gp],
                 qx_numerical_flux[gp],
                 qy_numerical_flux[gp]);
    }
}

void Land::GetEX(const Stepper& stepper,
                 const uint gp,
                 const Array2D<double>& surface_normal,
                 const std::vector<double>& ze_in,
                 const std::vector<double>& qx_in,
                 const std::vector<double>& qy_in,
                 double& ze_ex,
                 double& qx_ex,
                 double& qy_ex) {
    double n_x, n_y, t_x, t_y, qn_in, qt_in, qn_ex, qt_ex;

    n_x = surface_normal[gp][GlobalCoord::x];
    n_y = surface_normal[gp][GlobalCoord::y];
    t_x = -n_y;
    t_y = n_x;

    qn_in = qx_in[gp] * n_x + qy_in[gp] * n_y;
    qt_in = qx_in[gp] * t_x + qy_in[gp] * t_y;

    qn_ex = -qn_in;
    qt_ex = qt_in;

    ze_ex = ze_in[gp];
    qx_ex = qn_ex * n_x + qt_ex * t_x;
    qy_ex = qn_ex * n_y + qt_ex * t_y;
}
}
}

#endif