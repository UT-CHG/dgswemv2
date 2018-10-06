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
    HybMatrix<double, SWE::n_variables> q_ex;

  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    void Initialize(uint ngp);

    void ComputeFlux(const RKStepper& stepper,
                     const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                     const HybMatrix<double, SWE::n_variables>& q_in,
                     const HybMatrix<double, SWE::n_variables>& aux_in,
                     HybMatrix<double, SWE::n_variables>& F_hat);

    void ComputeFlux(const RKStepper& stepper,
                     const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                     const Column<HybMatrix<double, SWE::n_variables>>& q_in,
                     const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
                     Column<HybMatrix<double, SWE::n_variables>>&& F_hat);

    StatVector<double, SWE::n_variables> GetEX(const RKStepper& stepper,
                                               const StatVector<double, SWE::n_dimensions>& surface_normal,
                                               const StatVector<double, SWE::n_variables>& q_in);
};

template <typename BoundaryType>
void Land::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);
    this->q_ex.resize(SWE::n_variables, ngp);
}

void Land::Initialize(uint ngp) {
    this->q_ex.resize(SWE::n_variables, ngp);
}

void Land::ComputeFlux(const RKStepper& stepper,
                       const HybMatrix<double, SWE::n_dimensions>& surface_normal,
                       const HybMatrix<double, SWE::n_variables>& q_in,
                       const HybMatrix<double, SWE::n_variables>& aux_in,
                       HybMatrix<double, SWE::n_variables>& F_hat) {
    // *** //
    auto n_x = row(surface_normal, GlobalCoord::x);
    auto n_y = row(surface_normal, GlobalCoord::y);
    auto t_x = -n_y;
    auto t_y = n_x;

    auto qn_ex = -(vec_cw_mult(row(q_in, SWE::Variables::qx), n_x) + vec_cw_mult(row(q_in, SWE::Variables::qy), n_y));
    auto qt_ex = vec_cw_mult(row(q_in, SWE::Variables::qx), t_x) + vec_cw_mult(row(q_in, SWE::Variables::qy), t_y);

    row(this->q_ex, SWE::Variables::ze) = row(q_in, SWE::Variables::ze);
    row(this->q_ex, SWE::Variables::qx) = vec_cw_mult(qn_ex, n_x) + vec_cw_mult(qt_ex, t_x);
    row(this->q_ex, SWE::Variables::qy) = vec_cw_mult(qn_ex, n_y) + vec_cw_mult(qt_ex, t_y);

    for (uint gp = 0; gp < columns(q_in); ++gp) {
        column(F_hat, gp) = LLF_flux(
            Global::g, column(q_in, gp), column(this->q_ex, gp), column(aux_in, gp), column(surface_normal, gp));
    }
}

void Land::ComputeFlux(const RKStepper& stepper,
                       const Column<HybMatrix<double, SWE::n_dimensions>>& surface_normal,
                       const Column<HybMatrix<double, SWE::n_variables>>& q_in,
                       const Column<HybMatrix<double, SWE::n_auxiliaries>>& aux_in,
                       Column<HybMatrix<double, SWE::n_variables>>&& F_hat) {
    // *** //
    double n_x = surface_normal[GlobalCoord::x];
    double n_y = surface_normal[GlobalCoord::y];
    double t_x = -n_y;
    double t_y = n_x;

    double qn_ex = -(q_in[SWE::Variables::qx] * n_x + q_in[SWE::Variables::qy] * n_y);
    double qt_ex = q_in[SWE::Variables::qx] * t_x + q_in[SWE::Variables::qy] * t_y;

    this->q_ex(SWE::Variables::ze, 0) = q_in[SWE::Variables::ze];
    this->q_ex(SWE::Variables::qx, 0) = qn_ex * n_x + qt_ex * t_x;
    this->q_ex(SWE::Variables::qy, 0) = qn_ex * n_y + qt_ex * t_y;

    F_hat = LLF_flux(Global::g, q_in, column(this->q_ex, 0), aux_in, surface_normal);
}

StatVector<double, SWE::n_variables> Land::GetEX(const RKStepper& stepper,
                                                 const StatVector<double, SWE::n_dimensions>& surface_normal,
                                                 const StatVector<double, SWE::n_variables>& q_in) {
    StatVector<double, SWE::n_variables> q_ex;

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
