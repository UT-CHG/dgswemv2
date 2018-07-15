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
                     const std::vector<StatVector<double, SWE::n_variables>>& q_in,
                     const std::vector<StatVector<double, SWE::n_auxiliaries>>& aux_in,
                     std::vector<StatVector<double, SWE::n_variables>>& F_hat);

    void GetEX(const RKStepper& stepper,
               const std::vector<double>& surface_normal,
               const StatVector<double, SWE::n_variables>& q_in,
               StatVector<double, SWE::n_variables>& q_ex);
};

void Land::ComputeFlux(const RKStepper& stepper,
                       const Array2D<double>& surface_normal,
                       const std::vector<StatVector<double, SWE::n_variables>>& q_in,
                       const std::vector<StatVector<double, SWE::n_auxiliaries>>& aux_in,
                       std::vector<StatVector<double, SWE::n_variables>>& F_hat) {
    StatVector<double, SWE::n_variables> q_ex;
    for (uint gp = 0; gp < q_in.size(); ++gp) {
        this->GetEX(stepper, surface_normal[gp], q_in[gp], q_ex);

        LLF_flux(Global::g, q_in[gp], q_ex, aux_in[gp], surface_normal[gp], F_hat[gp]);
    }
}

void Land::GetEX(const RKStepper& stepper,
                 const std::vector<double>& surface_normal,
                 const StatVector<double, SWE::n_variables>& q_in,
                 StatVector<double, SWE::n_variables>& q_ex) {
    double n_x, n_y, t_x, t_y, qn_in, qt_in, qn_ex, qt_ex;

    n_x = surface_normal[GlobalCoord::x];
    n_y = surface_normal[GlobalCoord::y];
    t_x = -n_y;
    t_y = n_x;

    qn_in = q_in[SWE::Variables::qx] * n_x + q_in[SWE::Variables::qy] * n_y;
    qt_in = q_in[SWE::Variables::qx] * t_x + q_in[SWE::Variables::qy] * t_y;

    qn_ex = -qn_in;
    qt_ex = qt_in;

    q_ex[SWE::Variables::ze] = q_in[SWE::Variables::ze];
    q_ex[SWE::Variables::qx] = qn_ex * n_x + qt_ex * t_x;
    q_ex[SWE::Variables::qy] = qn_ex * n_y + qt_ex * t_y;
}
}
}
}

#endif