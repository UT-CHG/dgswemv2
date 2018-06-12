#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "../general_definitions.hpp"
#include "../preprocessor/input_parameters.hpp"

class Stepper {
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
    double ramp;

  public:
    Stepper() = default;
    Stepper(const StepperInput& stepper_input);

    uint GetOrder() const { return this->order; }
    uint GetNumStages() const { return this->nstages; }
    double GetDT() const { return this->dt; }

    uint GetStep() const { return this->step; }
    uint GetStage() const { return this->stage; }
    uint GetTimestamp() const { return this->timestamp; }

    double GetTimeAtCurrentStage() const { return this->t + this->dt * this->drk[this->stage]; }
    double GetRamp() const { return 1.0; /*tanh(this->GetTimeAtCurrentStage() / 86400.0);*/ }

    Stepper& operator++() {
        ++(this->stage);
        ++(this->timestamp);

        this->stage = this->stage % this->nstages;

        if (this->stage == 0) {
            this->t += this->dt;
            ++(this->step);
        }

        return *this;
    }
};

#endif