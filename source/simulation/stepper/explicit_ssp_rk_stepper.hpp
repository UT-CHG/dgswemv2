#ifndef E_SSP_RK_STEPPER_HPP
#define E_SSP_RK_STEPPER_HPP

#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "preprocessor/input_parameters.hpp"

/**
 * Explicit Strong Stability preserving Runge-Kutta methods
 * Class discretizes an ODE using strong stability preserving Runge Kutta methods.
 */
class ESSPRKStepper : public Stepper {
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
    ESSPRKStepper() = default;
    ESSPRKStepper(const StepperInput& stepper_input);

    uint GetOrder() const { return this->order; }
    uint GetNumStages() const { return this->nstages; }
    double GetDT() const { return this->dt; }

    void SetDT(double dt) { this->dt = dt; };

    uint GetStep() const { return this->step; }
    uint GetStage() const { return this->stage; }
    uint GetTimestamp() const { return this->timestamp; }

    double GetTimeAtCurrentStage() const { return this->t + this->dt * this->drk[this->stage]; }
    double GetRamp() const { return this->ramp; }

    ESSPRKStepper& operator++() {
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

        for ( uint var = 0; var < elt.data.state.size(); ++var ) {
            set_constant(next_state.q[var], 0.0);

            for (uint s = 0; s <= this->stage; ++s) {
                next_state.q[var] += this->ark[stage][s] * state[s].q[var] + this->dt * this->brk[stage][s] * state[s].solution[var];
            }

            if (this->stage + 1 == this->nstages)
                std::swap(state[0].q, state[this->nstages].q);
        }
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
void ESSPRKStepper::save(Archive& ar, unsigned) const {
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
void ESSPRKStepper::load(Archive& ar, unsigned) {
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