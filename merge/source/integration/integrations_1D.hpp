#ifndef INTEGRATIONS_1D_HPP
#define INTEGRATIONS_1D_HPP

#include "../general_definitions.hpp"

namespace Integration {
	class GaussLegendre_1D : Integration<1> {
	public:
		std::pair<std::vector<double>, std::vector<Point<1>>> get_rule(uint);
		void rule_test(uint, const std::pair<std::vector<double>, std::vector<Point<1>>>&);
	private:
		std::pair<std::vector<double>, std::vector<Point<1>>> gp_data(uint);
	};
}

#endif