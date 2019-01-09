#ifndef RKDG_SWE_BC_OUTFLOW_HPP
#define RKDG_SWE_BC_OUTFLOW_HPP

#include "problem/SWE/problem_function_files/swe_initial_condition_functions.hpp"

namespace SWE {
namespace RKDG {
namespace BC {
class Outflow {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound);

    template <typename StepperType, typename BoundaryType>
    void ComputeFlux(const StepperType& stepper, BoundaryType& bound);
};

template <typename BoundaryType>
void Outflow::Initialize(BoundaryType& bound) {}

template <typename StepperType, typename BoundaryType>
void Outflow::ComputeFlux(const StepperType& stepper, BoundaryType& bound) {
    auto& boundary = bound.data.boundary[bound.bound_id];

    for (uint gp = 0; gp < columns(boundary.q_at_gp); ++gp) {
        LLF_flux(Global::g,
                 column(boundary.q_at_gp, gp),
                 column(boundary.q_at_gp, gp),
                 column(boundary.aux_at_gp, gp),
                 column(bound.surface_normal, gp),
                 column(boundary.F_hat_at_gp, gp));
    }
}
}
}
}

#endif