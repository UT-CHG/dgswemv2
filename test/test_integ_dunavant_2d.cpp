#include "integration/integrations_2D.hpp"

int main() {
    bool any_error = false;
    Integration::Dunavant_2D dunavant;
    std::pair<std::vector<double>, std::vector<Point<2>>> rule;

    for (uint p = 1; p < 20; ++p) {
		    double exact_integration = 1 / ((double)p + 1)*((1 - pow(-1.0, p)) / ((double)p + 2) + 2 * pow(-1.0, p)); // S(x^p)dxdy over triangle

        rule = dunavant.GetRule(p);

		    double num_integration = 0;
    		for (uint gp = 0; gp < rule.first.size(); gp++) { 
            num_integration += pow(rule.second[gp][GlobalCoord::x], p)*rule.first[gp]; 
        }

        double err = abs((num_integration - exact_integration) / exact_integration);

        if (err > 100 * std::numeric_limits<double>::epsilon()) {
            any_error = true;
            std::cerr << 
                "Error found in Dunavant 2D at " << std::to_string(p) << 
            std::endl;
        }
    }

    if (any_error) {
        return 1;
    }

    return 0;
}