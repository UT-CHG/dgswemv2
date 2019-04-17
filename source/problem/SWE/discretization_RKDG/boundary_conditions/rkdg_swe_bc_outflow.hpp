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

    auto& n_x = bound.surface_normal[GlobalCoord::x];
    auto& n_y = bound.surface_normal[GlobalCoord::y];

    std::array<double,SWE::n_variables> F_hat_tmp;

    for (uint gp = 0; gp < boundary.q_at_gp[0].size(); ++gp) {
        LLF_flux(Global::g,
                                  boundary.q_at_gp[SWE::Variables::ze][gp],
                 boundary.q_at_gp[SWE::Variables::qx][gp],
                 boundary.q_at_gp[SWE::Variables::qy][gp],
                 boundary.q_at_gp[SWE::Variables::ze][gp],
                 boundary.q_at_gp[SWE::Variables::qx][gp],
                 boundary.q_at_gp[SWE::Variables::qy][gp],
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