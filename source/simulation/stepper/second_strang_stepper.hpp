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
        : first(stepper_input), second(stepper_input), current(&first) {
        this->first.SetDT(this->first.GetDT() / 2.0);
    }

    First& GetFirstStepper() { return this->first; }
    Second& GetSecondStepper() { return this->second; }

    SecondStrangStepper<First, Second>& operator=(SecondStrangStepper<First, Second>&& rhs) {
        this->first  = std::move(rhs.first);
        this->second = std::move(rhs.second);

        return *this;
    }
};

#endif