#include "simulation/stepper/rk_stepper.hpp"
#include "preprocessor/input_parameters.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

// This test checks the accuracy of the RKSSP methods by solving
// y''+y = t/2, whose solution is sin(t) + t/2;

int main() {
    bool error_found = false;

    std::array<std::pair<int, int>, 11> rk_pairs;
    rk_pairs[0]  = {1, 1};
    rk_pairs[1]  = {2, 2};
    rk_pairs[2]  = {3, 3};
    rk_pairs[3]  = {4, 3};
    rk_pairs[4]  = {5, 3};
    rk_pairs[5]  = {6, 3};
    rk_pairs[6]  = {7, 3};
    rk_pairs[7]  = {5, 4};
    rk_pairs[8]  = {6, 4};
    rk_pairs[9]  = {7, 4};
    rk_pairs[10] = {8, 4};

    double dt   = 0.00005;
    uint nsteps = 5. / dt + 1;

    using State      = std::array<double, 2>;
    auto compute_rhs = [](State y, double t) -> State { return {y[1], -y[0] + 0.5 * t}; };

    for (auto& pair : rk_pairs) {
        StepperInput stepper_input;

        stepper_input.nstages = pair.first;
        stepper_input.order   = pair.second;
        stepper_input.dt      = dt;

        RKStepper rk_stepper(stepper_input);

        double t = 0;

        std::vector<State> y(pair.first + 1);
        y[0] = {0, 1.5};
        std::vector<State> rhs(pair.first);

        for (uint step = 0; step < nsteps; ++step) {
            for (uint stage = 0; stage < rk_stepper.get_num_stages(); ++stage) {
                rhs[stage] = compute_rhs(y[stage], rk_stepper.get_t_at_curr_stage());
                y[stage + 1] = {0, 0};
                for (uint s = 0; s < stage + 1; ++s) {
                    y[stage + 1][0] += rk_stepper.ark[stage][s] * y[s][0] + dt * rk_stepper.brk[stage][s] * rhs[s][0];

                    y[stage + 1][1] += rk_stepper.ark[stage][s] * y[s][1] + dt * rk_stepper.brk[stage][s] * rhs[s][1];
                }
                ++rk_stepper;
            }
            std::swap(y[0], y[rk_stepper.get_num_stages()]);
        }

        double t = rk_stepper.get_t_at_curr_stage();

        std::cout << "At time: " << t << "\n";
        std::cout << "Got: " << std::setprecision(14) << y[0][0] << " Should be: " << std::sin(t) + 0.5 * t << "\n";
        std::cout << "Got: " << y[0][1] << " Should be: " << std::cos(t) + 0.5 << "\n\n";

        if (std::abs(y[0][0] - std::sin(t) - 0.5 * t) > std::pow(10., -2)) {
            std::cerr << "Error in Runge-Kutta timestepping routine\n";
            std::cerr << "RK scheme does not seem to reproduce accurate results\n";
            error_found = true;
        }
    }

    if (error_found) {
        return 1;
    }
    return 0;
}
