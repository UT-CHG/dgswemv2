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

    double t = stepper.GetTimeAtCurrentStage();

    this->q_ex = bound.ComputeFgp([t](Point<2>& pt) {
        double ze = 0.0;
        double qx = 0.0;
        double qy = 0.0;

        if (t <= 3.0) {
            ze = cos(PI * t) - 1.0;
        } else {
            ze = -2.0;
        }

        StatVector<double, SWE::n_variables> q{ze, qx, qy};
        // StatVector<double, SWE::n_variables> q(SWE::ic_q(t, pt));

        return q;
    });

    for (uint gp = 0; gp < columns(boundary.q_at_gp); ++gp) {
        LLF_flux(Global::g,
                 column(boundary.q_at_gp, gp),
                 column(this->q_ex, gp),
                 column(boundary.aux_at_gp, gp),
                 column(bound.surface_normal, gp),
                 column(boundary.F_hat_at_gp, gp));
    }
}
}
}
}

#endif