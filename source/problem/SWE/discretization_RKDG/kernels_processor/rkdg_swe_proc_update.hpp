#ifndef RKDG_SWE_PROC_UPDATE_HPP
#define RKDG_SWE_PROC_UPDATE_HPP

#include "problem/SWE/discretization_RKDG/rkdg_swe_problem.hpp"

namespace SWE {
namespace RKDG {

class UpdateKernel {
private:
    using StepperType = typename Problem::ProblemStepperType;
public:
    constexpr static bool is_vectorized = true;

    UpdateKernel(StepperType& stepper_) : stepper(stepper_) {}

    template <typename SoA>
    void operator() (SoA& soa) const {
        auto& state = soa.data.state[stepper.GetStage()];

        state.solution[SWE::Variables::ze] = soa.ApplyMinv(state.rhs[SWE::Variables::ze]);
        state.solution[SWE::Variables::qx] = soa.ApplyMinv(state.rhs[SWE::Variables::qx]);
        state.solution[SWE::Variables::qy] = soa.ApplyMinv(state.rhs[SWE::Variables::qy]);

        stepper.UpdateState(soa);
    }


private:
    StepperType& stepper;
};

}
}

#endif