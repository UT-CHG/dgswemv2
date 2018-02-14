#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "general_definitions.hpp"

struct Stepper {
    Stepper()=default;
    Stepper(uint nstages, uint order, double dt, double T_end, double t0 = 0.);

    void InitializeCoefficients();

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

    uint order;
    uint nstages;
    uint irk;

    std::vector<std::vector<double>> ark;
    std::vector<std::vector<double>> brk;
    std::vector<std::vector<double>> crk;
    std::vector<double> drk;

    uint step;
    uint timestamp;
    double _t;
    double _dt;
    double T_end;

#ifdef HAS_HPX
    template <typename Archive>
    void save(Archive& ar, unsigned) const;

    template <typename Archive>
    void load(Archive& ar, unsigned);

    HPX_SERIALIZATION_SPLIT_MEMBER()
#endif
};

#ifdef HAS_HPX
template<typename Archive>
void Stepper::save(Archive& ar, unsigned) const {
    ar & order & nstages & irk & _dt & timestamp & _t & T_end;
}

template<typename Archive>
void Stepper::load(Archive& ar, unsigned) {
    ar & order & nstages & irk & _dt & timestamp & _t & T_end;

    step = timestamp/nstages;
    InitializeCoefficients();
}
#endif
#endif