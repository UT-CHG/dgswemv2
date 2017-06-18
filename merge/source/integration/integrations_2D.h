#ifndef INTEGRATIONS_2D_H
#define INTEGRATIONS_2D_H

#include "../general_definitions.h"

namespace Integration{
class Dunavant_2D : Integration<2> {
public:
	std::pair<std::vector<double>, std::vector<Point<2>>> get_rule(int);
	void rule_test(int, const std::pair<std::vector<double>, std::vector<Point<2>>>&);
private:
	std::vector<int> permutation_data(int);
	std::pair<std::vector<double>, std::vector<Point<3>>> gp_data(int);
};
}

#endif