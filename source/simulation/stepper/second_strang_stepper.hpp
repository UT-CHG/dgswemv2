#ifndef SECOND_STRANG_STEPPER_HPP
#define SECOND_STRANG_STEPPER_HPP

#include "general_definitions.hpp"
#include "preprocessor/input_parameters.hpp"

/**
 * Second order strang splitting stepper
 */
template <typename First, typename Second>
class SecondStrangStepper {
  private:
    First first;
    Second second;

  public:
    SecondStrangStepper() = default;
    SecondStrangStepper(const StepperInput& stepper_input)
        : first(stepper_input), second(stepper_input) {
        this->first.SetDT(this->first.GetDT() / 2.0);
    }

    First& GetFirstStepper() { return this->first; }
    Second& GetSecondStepper() { return this->second; }

    const First& GetFirstStepper() const { return this->first; }
    const Second& GetSecondStepper() const { return this->second; }

    // In second strang spliting it's the second stepper that does major stepping
    uint GetStep() const { return this->second.GetStep(); }
    double GetTimeAtCurrentStage() const { return this->second.GetTimeAtCurrentStage(); }

    SecondStrangStepper<First, Second>& operator=(SecondStrangStepper<First, Second>&& rhs) {
        this->first  = std::move(rhs.first);
        this->second = std::move(rhs.second);

        return *this;
    }
};

#endif