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

        Utilities::ignore(ze, qx, qy);

        if (t <= 3.0) {
            ze = cos(PI * t) - 1.0;
        } else {
            ze = -2.0;
        }

        // StatVector<double, SWE::n_variables> q{ze, qx, qy};
        StatVector<double, SWE::n_variables> q(SWE::ic_q(t, pt));

        return q;
    });

    for (uint gp = 0; gp < bound.data.get_ngp_boundary(bound.bound_id); ++gp) {
        NCP_flux(Global::g,
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
        const double un =
            roe_un(bound.surface_normal(GlobalCoord::x, gp),
                   bound.surface_normal(GlobalCoord::y, gp),
                   boundary.q_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp),
                   boundary.q_at_gp(SWE::Variables::qx, gp),
                   boundary.q_at_gp(SWE::Variables::qy, gp),
                   this->q_ex(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp),
                   this->q_ex(SWE::Variables::qx, gp),
                   this->q_ex(SWE::Variables::qy, gp));
        if (Utilities::almost_equal(un, 0.0)) {
            boundary.qb_hat_at_gp[gp] = 0.0;
        } else if (un > 0.0) {
            boundary.qb_hat_at_gp[gp] =
                transpose(column(bound.surface_normal, gp)) *
                bed_flux(boundary.q_at_gp(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp),
                         boundary.q_at_gp(SWE::Variables::qx, gp),
                         boundary.q_at_gp(SWE::Variables::qy, gp));
        } else if (un < 0.0) {
            boundary.qb_hat_at_gp[gp] =
                transpose(column(bound.surface_normal, gp)) *
                bed_flux(this->q_ex(SWE::Variables::ze, gp) + boundary.aux_at_gp(SWE::Auxiliaries::bath, gp),
                         this->q_ex(SWE::Variables::qx, gp),
                         this->q_ex(SWE::Variables::qy, gp));
        }
    }
}
}
}
}

#endif