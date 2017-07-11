#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "../general_definitions.hpp"

struct Stepper {
	Stepper(uint nstages, uint order, double dt);

	inline
		uint get_num_stages() const
	{
		return nstages;
	}

	inline
		double get_t_at_curr_stage() const
	{
		return _t + _dt*drk[irk];
	}

	inline
		double get_dt() const
	{
		return _dt;
	}

	inline
		uint get_stage() const
	{
		return irk;
	}

	inline
		Stepper& operator++()
	{
		++irk;
		irk = irk % nstages;

		if (irk == 0) {
			_t += _dt;
		}

		return *this;
	}

	const uint nstages;
	uint irk;

	std::vector<std::vector<double> > ark;
	std::vector<std::vector<double> > brk;
	std::vector<std::vector<double> > crk;
	std::vector<double> drk;

	double _t;
	const double _dt;
};

#endif
