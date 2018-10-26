#ifndef RK_STEPPER_HPP
#define RK_STEPPER_HPP

#include "utilities/almost_equal.hpp"
#include "general_definitions.hpp"
#include "preprocessor/input_parameters.hpp"

/**
 * Strong Stability preserving Runge-Kutta methods
 * Class discretizes an ODE using strong stability preserving Runge Kutta methods.
 */
class RKStepper {
  public:
    Array2D<double> ark;
    Array2D<double> brk;
    Array2D<double> crk;
    std::vector<double> drk;

  private:
    uint order;
    uint nstages;
    double dt;

    uint step;
    uint stage;
    uint timestamp;

    double t;
    double ramp_duration;
    double ramp;

  public:
    RKStepper() = default;
    RKStepper(const StepperInput& stepper_input);

    /**
     * Get the order accuracy of the RK stepper
     */
    uint GetOrder() const { return this->order; }

    /**
     * Get the number of stages per timestep
     */
    uint GetNumStages() const { return this->nstages; }

    /**
     * Get the size of the timestep (in seconds)
     */
    double GetDT() const { return this->dt; }

    /**
     * Get the current step number
     */
    uint GetStep() const { return this->step; }

    /**
     * Get the current stage number
     */
    uint GetStage() const { return this->stage; }

    /**
     * Get the current timestamp
     * The timestamp is the total number of stages that have been executed up to this point.
     */
    uint GetTimestamp() const { return this->timestamp; }

    /**
     * Get the simulated time at the current step and stage
     */
    double GetTimeAtCurrentStage() const { return this->t + this->dt * this->drk[this->stage]; }

    /**
     * Get ramp factor
     * Often for stability reasons, we scale boundary or source terms by a number that goes from 0 to 1
     * as the simulation begins. We refer to this factor as ramp. In this stepper, we
     * use a hyperbolic tangent function. When `t = ramp_duration`, `ramp = tanh(2)`.
     */
    double GetRamp() const { return this->ramp; }

    /**
     * Prefix incrementor advances the stepper by one RK-stage
     * This operation will update all of the internal states of the stepper by one Runge-Kutta stage.
     */
    RKStepper& operator++() {
        ++(this->stage);
        ++(this->timestamp);

        this->stage = this->stage % this->nstages;

        if (this->stage == 0) {
            this->t += this->dt;
            ++(this->step);
        }

        if (!Utilities::almost_equal(this->ramp_duration, 0)) {
            this->ramp = std::tanh(2 * (this->GetTimeAtCurrentStage() / 86400) / this->ramp_duration);
        }

        return *this;
    }

    template <typename ElementType>
    void UpdateState(ElementType& elt) const {
        auto& state      = elt.data.state;
        auto& next_state = elt.data.state[this->stage + 1];

        set_constant(next_state.q, 0.0);

        for (uint s = 0; s <= this->stage; ++s) {
            next_state.q += this->ark[stage][s] * state[s].q + this->dt * this->brk[stage][s] * state[s].solution;
        }

        if (this->stage + 1 == this->nstages)
            std::swap(state[0].q, state[this->nstages].q);
    }

#ifdef HAS_HPX
    template <typename Archive>
    void save(Archive& ar, unsigned) const;

    template <typename Archive>
    void load(Archive& ar, unsigned);

    HPX_SERIALIZATION_SPLIT_MEMBER()
#endif
  private:
    void InitializeCoefficients();
};
#ifdef HAS_HPX
template <typename Archive>
void RKStepper::save(Archive& ar, unsigned) const {
    // clang-format off
    ar & order
       & nstages
       & dt
       & stage
       & timestamp
       & t
       & ramp_duration;
    // clang-format on
}

template <typename Archive>
void RKStepper::load(Archive& ar, unsigned) {
    // clang-format off
    ar & order
       & nstages
       & dt
       & stage
       & timestamp
       & t
       & ramp_duration;
    // clang-format on

    step = timestamp / nstages;
    InitializeCoefficients();

    if (!Utilities::almost_equal(this->ramp_duration, 0)) {
        this->ramp = std::tanh(2 * (this->GetTimeAtCurrentStage() / 86400) / this->ramp_duration);
    } else {
        this->ramp = 1;
    }
}
#endif
#endif