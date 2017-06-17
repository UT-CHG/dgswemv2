#ifndef INTEGRATION_RULES_1D_H
#define INTEGRATION_RULES_1D_H

#include "../general_definitions.h"

namespace Integration{
class GaussLegendre_1D : Integration<1> {
public:
	std::pair<std::vector<double>, std::vector<Point<1>>> get_rule(int);
	void rule_test(int, const std::pair<std::vector<double>, std::vector<Point<1>>>&);
private:
	std::pair<std::vector<double>, std::vector<Point<1>>> gp_data(int);
};
}

#endif