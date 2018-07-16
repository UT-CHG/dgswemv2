#include "integration/integrations_1D.hpp"
#include "utilities/almost_equal.hpp"

int main() {
    using Utilities::almost_equal;

    bool any_error = false;
    Integration::GaussLegendre_1D gausslegendre;
    std::pair<DynVector<double>, DynVector<Point<1>>> rule;

    for (uint p = 1; p < 20; ++p) {
        double exact_integration = 2 - 1 / ((double)p + 1) * (1 + pow(-1.0, p));  // S(1-x^p)dx from -1 to 1

        rule = gausslegendre.GetRule(p);

        double num_integration = 0;
        for (uint gp = 0; gp < rule.first.size(); gp++) {
            num_integration += (1.0 - pow(rule.second[gp][GlobalCoord::x], p)) * rule.first[gp];
        }

        if (!almost_equal(num_integration, exact_integration, 1.e+03)) {
            any_error = true;
            std::cerr << "Error found in Gauss-Legendre 1D at " << std::to_string(p)
                      << " - integration true value: " << exact_integration
                      << ", integration computed value: " << num_integration << std::endl;
        }
    }

    if (any_error) {
        return 1;
    }

    return 0;
}