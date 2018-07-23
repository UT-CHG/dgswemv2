#ifndef RKDG_SWE_BC_LAND_HPP
#define RKDG_SWE_BC_LAND_HPP

#include "general_definitions.hpp"
#include "simulation/stepper/rk_stepper.hpp"
#include "problem/SWE/discretization_RKDG/numerical_fluxes/rkdg_swe_numerical_fluxes.hpp"

namespace SWE {
namespace RKDG {
namespace BC {
class Land {
  private:
    DynMatrix<double> q_ex;

  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void ComputeFlux(const RKStepper& stepper,
                     const DynMatrix<double>& surface_normal,
                     const DynMatrix<double>& q_in,
                     const DynMatrix<double>& aux_in,
                     DynMatrix<double>& F_hat);

    DynVector<double> GetEX(const RKStepper& stepper,
                            const DynVector<double>& surface_normal,
                            const DynVector<double>& q_in);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);
    this->q_ex.resize(SWE::n_variables, ngp);
}

void Land::ComputeFlux(const RKStepper& stepper,
                       const DynMatrix<double>& surface_normal,
                       const DynMatrix<double>& q_in,
                       const DynMatrix<double>& aux_in,
                       DynMatrix<double>& F_hat) {
    // *** //
    auto n_x = row(surface_normal, GlobalCoord::x);
    auto n_y = row(surface_normal, GlobalCoord::y);
    auto t_x = -n_y;
    auto t_y = n_x;

    auto qn_ex = -(cwise_multiplication(row(q_in, SWE::Variables::qx), n_x) +
                   cwise_multiplication(row(q_in, SWE::Variables::qy), n_y));
    auto qt_ex = cwise_multiplication(row(q_in, SWE::Variables::qx), t_x) +
                 cwise_multiplication(row(q_in, SWE::Variables::qy), t_y);

    row(this->q_ex, SWE::Variables::ze) = row(q_in, SWE::Variables::ze);
    row(this->q_ex, SWE::Variables::qx) = cwise_multiplication(qn_ex, n_x) + cwise_multiplication(qt_ex, t_x);
    row(this->q_ex, SWE::Variables::qy) = cwise_multiplication(qn_ex, n_y) + cwise_multiplication(qt_ex, t_y);

    for (uint gp = 0; gp < columns(q_in); ++gp) {
        column(F_hat, gp) = LLF_flux(
            Global::g, column(q_in, gp), column(this->q_ex, gp), column(aux_in, gp), column(surface_normal, gp));
    }
}

DynVector<double> Land::GetEX(const RKStepper& stepper,
                              const DynVector<double>& surface_normal,
                              const DynVector<double>& q_in) {
    DynVector<double> q_ex(SWE::n_variables);

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

    return q_ex;
}
}
}
}

#endif