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
                     const DynVector<StatVector<double, SWE::n_dimensions>>& surface_normal,
                     const DynMatrix<double>& q_in,
                     const DynMatrix<double>& aux_in,
                     DynMatrix<double>& F_hat);

    void GetEX(const RKStepper& stepper,
               const StatVector<double, SWE::n_dimensions>& surface_normal,
               const DynVector<double>& q_in,
               DynVector<double>& q_ex);
};

void Land::ComputeFlux(const RKStepper& stepper,
                       const DynVector<StatVector<double, SWE::n_dimensions>>& surface_normal,
                       const DynMatrix<double>& q_in,
                       const DynMatrix<double>& aux_in,
                       DynMatrix<double>& F_hat) {
    // *** //
    DynVector<double> q_ex(SWE::n_variables);
    for (uint gp = 0; gp < columns(q_in); ++gp) {
        this->GetEX(stepper, surface_normal[gp], column(q_in, gp), q_ex);

        column(F_hat, gp) = LLF_flux(Global::g, column(q_in, gp), q_ex, column(aux_in, gp), surface_normal[gp]);
    }
}

void Land::GetEX(const RKStepper& stepper,
                 const StatVector<double, SWE::n_dimensions>& surface_normal,
                 const DynVector<double>& q_in,
                 DynVector<double>& q_ex) {
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