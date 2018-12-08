#ifndef RKDG_SWE_BC_FUNCTION_HPP
#define RKDG_SWE_BC_FUNCTION_HPP

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"

namespace SWE {
namespace RKDG {
namespace BC {
class Function {
  private:
    HybMatrix<double, SWE::n_variables> q_ex;

  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename BoundaryType>
    void ComputeFlux(const StepperType& stepper, BoundaryType& bound);
};

template <typename BoundaryType>
void Function::Initialize(BoundaryType& bound) {
    uint ngp = bound.data.get_ngp_boundary(bound.bound_id);
    this->q_ex.resize(SWE::n_variables, ngp);
}

template <typename StepperType, typename BoundaryType>
void Function::ComputeFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];

    uint ngp = boundary.F_hat_at_gp[0].size();
    double t    = stepper.GetTimeAtCurrentStage();
    auto q_func = [t](Point<2>& pt) { return SWE::ic_q(t, pt); };

    this->q_ex = bound.ComputeFgp(q_func);

    auto n_x = row(bound.surface_normal, GlobalCoord::x);
    auto n_y = row(bound.surface_normal, GlobalCoord::y);

    for (uint gp = 0; gp < ngp; ++gp) {

        std::array<double,SWE::n_variables> F_hat_tmp;

        LLF_flux(Global::g,
                 boundary.q_at_gp[SWE::Variables::ze][gp],
                 boundary.q_at_gp[SWE::Variables::qx][gp],
                 boundary.q_at_gp[SWE::Variables::qy][gp],
                 this->q_ex(SWE::Variables::ze,gp),
                 this->q_ex(SWE::Variables::qx,gp),
                 this->q_ex(SWE::Variables::qy,gp),
                 std::array<double,SWE::n_auxiliaries>{boundary.aux_at_gp[SWE::Auxiliaries::bath][gp],
                         boundary.aux_at_gp[SWE::Auxiliaries::h][gp],
                         boundary.aux_at_gp[SWE::Auxiliaries::sp][gp]},
                 std::array<double,SWE::n_dimensions>{n_x[gp],n_y[gp]},
                 F_hat_tmp
            );
        for ( uint var = 0; var < SWE::n_variables; ++var ) {
            boundary.F_hat_at_gp[var][gp] = F_hat_tmp[var];
        }
    }
}
}
}
}

#endif