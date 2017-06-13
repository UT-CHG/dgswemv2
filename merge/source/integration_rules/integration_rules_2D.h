#ifndef INTEGRATION_RULES_2D_H
#define INTEGRATION_RULES_2D_H

#include "../general_definitions.h"

class Dunavant_2D{
public:
	std::pair<std::vector<double>, std::vector<Point<2>>> get_rule(int);
	void rule_test(int, const std::pair<std::vector<double>, std::vector<Point<2>>>&);
private:
	std::vector<int> permutation_data(int);
	std::pair<std::vector<double>, std::vector<Point<3>>> gp_data(int);
};

#endif