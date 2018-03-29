#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "../general_definitions.hpp"

class Stepper {
  public:
    Array2D<double>     ark;
    Array2D<double>     brk;
    Array2D<double>     crk;
    std::vector<double> drk;

  private:
    const uint   nstages;
    const double dt;

    uint step;
    uint stage;
    uint timestamp;

    double t;

  public:
    Stepper(uint nstages, uint order, double dt);

    uint   GetNumStages() const { return this->nstages; }
    double GetDT() const { return this->dt; }

    uint   GetStep() const { return this->step; }
    uint   GetTimestamp() const { return this->timestamp; }
    uint   GetStage() const { return this->stage; }
    double GetTimeAtCurrentStage() const { return this->t + this->dt * this->drk[this->stage]; }

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