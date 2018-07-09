#ifndef RKDG_SWE_BC_LAND_HPP
#define RKDG_SWE_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {} /*nothing to initialize*/

    void ComputeFlux(const RKStepper& stepper,
                     const Array2D<double>& surface_normal,
                     const std::vector<double>& sp_in,
                     const std::vector<double>& bath_in,
                     const std::vector<Vector<double, SWE::n_variables>>& u_in,
                     std::vector<Vector<double, SWE::n_variables>>& F_hat);

    void GetEX(const RKStepper& stepper,
               const std::vector<double>& surface_normal,
               const Vector<double, SWE::n_variables>& u_in,
               Vector<double, SWE::n_variables>& u_ex);
};

void Land::ComputeFlux(const RKStepper& stepper,
                       const Array2D<double>& surface_normal,
                       const std::vector<double>& sp_in,
                       const std::vector<double>& bath_in,
                       const std::vector<Vector<double, SWE::n_variables>>& u_in,
                       std::vector<Vector<double, SWE::n_variables>>& F_hat) {
    Vector<double, SWE::n_variables> u_ex;
    for (uint gp = 0; gp < u_in.size(); ++gp) {
        this->GetEX(stepper, surface_normal[gp], u_in[gp], u_ex);

        LLF_flux(Global::g, u_in[gp], u_ex, bath_in[gp], sp_in[gp], surface_normal[gp], F_hat[gp]);
    }
}

void Land::GetEX(const RKStepper& stepper,
                 const std::vector<double>& surface_normal,
                 const Vector<double, SWE::n_variables>& u_in,
                 Vector<double, SWE::n_variables>& u_ex) {
    double n_x, n_y, t_x, t_y, qn_in, qt_in, qn_ex, qt_ex;

    n_x = surface_normal[GlobalCoord::x];
    n_y = surface_normal[GlobalCoord::y];
    t_x = -n_y;
    t_y = n_x;

    qn_in = u_in[SWE::Variables::qx] * n_x + u_in[SWE::Variables::qy] * n_y;
    qt_in = u_in[SWE::Variables::qx] * t_x + u_in[SWE::Variables::qy] * t_y;

    qn_ex = -qn_in;
    qt_ex = qt_in;

    u_ex[SWE::Variables::ze] = u_in[SWE::Variables::ze];
    u_ex[SWE::Variables::qx] = qn_ex * n_x + qt_ex * t_x;
    u_ex[SWE::Variables::qy] = qn_ex * n_y + qt_ex * t_y;
}
}
}
}

#endif