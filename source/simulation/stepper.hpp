#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "../general_definitions.hpp"

class Stepper {
  public:
    Array2D<double> ark;
    Array2D<double> brk;
    Array2D<double> crk;
    std::vector<double> drk;

  private:
    const uint nstages;
    const double dt;

    uint step;
    uint stage;
    uint timestamp;

    double t;

  public:
    Stepper(uint nstages, uint order, double dt);

    uint GetNumStages() const { return nstages; }
    double GetDT() const { return dt; }
    
    uint GetStep() const { return step; }
    uint GetTimestamp() const { return timestamp; }
    uint GetStage() const { return stage; }
    double GetTimeAtCurrentStage() const { return t + dt * drk[stage]; }

    Stepper& operator++() {
        ++stage;
        ++timestamp;

        stage = stage % nstages;

        if (stage == 0) {
            t += dt;
            ++step;
        }

        return *this;
    }
};

#endif