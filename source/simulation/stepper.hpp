#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "general_definitions.hpp"

struct Stepper {
    Stepper(uint nstages, uint order, double dt);

    uint get_num_stages() const { return nstages; }

    double get_t_at_curr_stage() const { return _t + _dt * drk[irk]; }

    double get_dt() const { return _dt; }

    uint get_step() const { return step; }

    uint get_timestamp() const { return timestamp; }

    uint get_stage() const { return irk; }

    Stepper& operator++() {
        ++irk;
        ++timestamp;

        irk = irk % nstages;

        if (irk == 0) {
            _t += _dt;
            ++step;
        }

        return *this;
    }

    const uint nstages;
    uint irk;

    std::vector<std::vector<double>> ark;
    std::vector<std::vector<double>> brk;
    std::vector<std::vector<double>> crk;
    std::vector<double> drk;

    uint step;
    uint timestamp;
    double _t;
    const double _dt;
};

#endif