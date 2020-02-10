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

    template <typename StepperType, typename BoundaryType>
    void ComputeBedFlux(const StepperType& stepper, BoundaryType& bound);
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
        double hc = 0.0;
        Utilities::ignore(ze, qx, qy, hc);
        // return StatVector<double, SWE::n_variables>{ze, qx, qy, hc};
        return StatVector<double, SWE::n_variables>(SWE::ic_q(t, pt));
    });

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        HLL_flux(Global::g,
                 column(boundary.q_at_gp, gp),
                 column(this->q_ex, gp),
                 column(boundary.aux_at_gp, gp),
                 column(boundary.aux_at_gp, gp),
                 column(bound.surface_normal, gp),
                 column(boundary.F_hat_at_gp, gp));
    }
}

template <typename StepperType, typename BoundaryType>
void Function::ComputeBedFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        const double un = roe_un(column(boundary.q_at_gp, gp),
                                 column(this->q_ex, gp),
                                 column(boundary.aux_at_gp, gp),
                                 column(boundary.aux_at_gp, gp),
                                 column(bound.surface_normal, gp));
        if (Utilities::almost_equal(un, 0.0)) {
            boundary.qb_hat_at_gp[gp] = 0.0;
        } else if (un > 0.0) {
            boundary.qb_hat_at_gp[gp] = transpose(column(bound.surface_normal, gp)) *
                                        bed_flux(column(boundary.q_at_gp, gp), column(boundary.aux_at_gp, gp));
        } else if (un < 0.0) {
            boundary.qb_hat_at_gp[gp] = transpose(column(bound.surface_normal, gp)) *
                                        bed_flux(column(this->q_ex, gp), column(boundary.aux_at_gp, gp));
        }
    }
}
}
}
}

#endif