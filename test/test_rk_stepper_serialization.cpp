#include "utilities/almost_equal.hpp"
#include "simulation/stepper/rk_stepper.hpp"

// This test checks the accuracy of the RKSSP methods by solving
// y''+y = t/2, whose solution is sin(t) + t/2;

using State                = std::array<double, 2>;
constexpr auto compute_rhs = [](State y, double t) -> State { return {y[1], -y[0] + 0.5 * t}; };

State solve_ode(RKStepper& rk_stepper, const State& y0, uint nsteps) {
    std::vector<State> y(rk_stepper.GetNumStages() + 1);
    y[0] = y0;
    std::vector<State> rhs(rk_stepper.GetNumStages());

    for (uint step = 0; step < nsteps; ++step) {
        for (uint stage = 0; stage < rk_stepper.GetNumStages(); ++stage) {
            rhs[stage]   = compute_rhs(y[stage], rk_stepper.GetTimeAtCurrentStage());
            y[stage + 1] = {0, 0};
            for (uint s = 0; s < stage + 1; ++s) {
                y[stage + 1][0] +=
                    rk_stepper.ark[stage][s] * y[s][0] + rk_stepper.GetDT() * rk_stepper.brk[stage][s] * rhs[s][0];

                y[stage + 1][1] +=
                    rk_stepper.ark[stage][s] * y[s][1] + rk_stepper.GetDT() * rk_stepper.brk[stage][s] * rhs[s][1];
            }
            ++rk_stepper;
        }
        std::swap(y[0], y[rk_stepper.GetNumStages()]);
    }

    return y[0];
}

int main() {
    bool error_found{false};

    State y0{0, 1.5};
    double dt   = 0.00005;
    uint nsteps = 5. / dt + 1;

    StepperInput stepper_input;
    stepper_input.nstages = 2;
    stepper_input.order   = 2;
    stepper_input.dt      = dt;

    RKStepper rk1(stepper_input);
    State y_non_stop = solve_ode(rk1, y0, nsteps);

    RKStepper o_rk2(stepper_input);
    State y_interim = solve_ode(o_rk2, y0, nsteps / 2);

    std::vector<char> buffer;
    hpx::serialization::output_archive o_archive(buffer);
    o_archive << o_rk2;

    hpx::serialization::input_archive i_archive(buffer);
    RKStepper i_rk2;
    i_archive >> i_rk2;

    State y_stop = solve_ode(i_rk2, y_interim, nsteps - (nsteps / 2));

    std::cout << "Got non-stop: " << std::setprecision(14) << y_non_stop[0] << " Serialized: " << y_stop[0] << "\n";
    std::cout << "Got non-stop: " << std::setprecision(14) << y_non_stop[1] << " Serialized: " << y_stop[1] << "\n";

    error_found =
        !(Utilities::almost_equal(y_non_stop[0], y_stop[0]) && Utilities::almost_equal(y_non_stop[1], y_stop[1]));

    if (error_found) {
        return 1;
    }
    return 0;
}