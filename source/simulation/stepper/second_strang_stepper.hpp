#ifndef SECOND_STRANG_STEPPER_HPP
#define SECOND_STRANG_STEPPER_HPP

#include "general_definitions.hpp"
#include "preprocessor/input_parameters.hpp"

/**
 * Second order strang splitting stepper
 */
template <typename First, typename Second>
class SecondStrangStepper : public Stepper {
  private:
    First first;
    Second second;

    bool done_first  = false;
    bool done_second = false;

    Stepper* current = nullptr;

  public:
    SecondStrangStepper() = default;
    SecondStrangStepper(const StepperInput& stepper_input)
        : first(stepper_input), second(stepper_input), current(&first) {
        this->first.SetDT(this->first.GetDT() / 2.0);
    }

    SecondStrangStepper<First, Second>& operator=(SecondStrangStepper<First, Second>&& rhs) {
        this->first  = std::move(rhs.first);
        this->second = std::move(rhs.second);

        this->done_first  = false;
        this->done_second = false;
        this->current     = &(this->first);

        return *this;
    }

    uint GetOrder() const { return this->current->GetOrder(); }
    uint GetNumStages() const { return this->current->GetNumStages(); }
    double GetDT() const { return this->current->GetDT(); }

    void SetDT(double dt) { this->current->SetDT(dt); }

    uint GetStep() const { return this->current->GetStep(); }
    uint GetStage() const { return this->current->GetStage(); }
    uint GetTimestamp() const { return this->current->GetTimestamp(); }

    double GetTimeAtCurrentStage() const { return this->current->GetTimeAtCurrentStage(); }
    double GetRamp() const { return this->current->GetRamp(); }

    SecondStrangStepper<First, Second>& operator++() {
        ++(*this->current);

        if (this->current->GetStage() == 0) {
            if (!done_first && !done_second) {
                this->current = &this->second;
                done_first    = true;
            } else if (done_first && !done_second) {
                this->current = &this->first;
                done_second   = true;
            } else if (done_first && done_second) {
                this->current = &this->first;
                done_first    = false;
                done_second   = false;
            }
        }

        return *this;
    }

    template <typename ElementType>
    void UpdateState(ElementType& elt) const {
        if (!done_first && !done_second) {
            this->first.UpdateState(elt);
        } else if (done_first && !done_second) {
            this->second.UpdateState(elt);
        } else if (done_first && done_second) {
            this->first.UpdateState(elt);
        }
    }
};

#endif